# Adaptive Layer Height Algorithm in PrusaSlicer

## Overview

PrusaSlicer implements an **adaptive layer height** feature that automatically adjusts the layer thickness
based on the local surface geometry of the 3D model. Surfaces with steep slopes get thinner layers
(for better surface quality), while surfaces with gentle slopes get thicker layers (for faster printing).

## Key Source Files

| File | Purpose |
|------|---------|
| `src/libslic3r/SlicingAdaptive.hpp` | Class declaration for the adaptive slicing core |
| `src/libslic3r/SlicingAdaptive.cpp` | Core algorithm implementation |
| `src/libslic3r/Slicing.hpp` | `SlicingParameters` struct and layer height profile API |
| `src/libslic3r/Slicing.cpp` | `layer_height_profile_adaptive()` – profile generation entry point |
| `src/slic3r/GUI/GLCanvas3D.cpp` | GUI integration – layer height editor with quality slider |

## Theoretical Background

The algorithm is based on the paper:

> **Florens Wasserfall, Norman Hendrich, Jianwei Zhang**  
> *Adaptive Slicing for the FDM Process Revisited*  
> 13th IEEE Conference on Automation Science and Engineering (CASE-2017),
> August 20–23, Xi'an, China. DOI: 10.1109/COASE.2017.8256074  
> https://tams.informatik.uni-hamburg.de/publications/2017/Adaptive%20Slicing%20for%20the%20FDM%20Process%20Revisited.pdf

The core idea is to **minimize surface roughness** (the "stairstepping" artifact) while maximizing
layer height to keep print time reasonable. Roughness is approximated as a geometric error triangle
whose area depends on the layer height **and** the local slope of the surface.

### Surface Roughness Model

When printing a surface that makes an angle `α` with the horizontal plane:

- **Horizontal surface (α = 0°)**: Layer stacking creates no visible stairstep artifacts.
- **45° slope**: Each layer creates a stairstep of height H × tan(α).
- **Vertical surface (α = 90°)**: Any layer height results in maximum stairstep artifacts.

The **error triangle** formed between the ideal surface and the actual stairstepped surface has area
proportional to `H² × sin(α) / cos(α)`, where `H` is the layer height. To keep roughness constant,
we solve for `H`:

```
H = C × sqrt(cos(α) / sin(α))
  = C × sqrt(n_cos / n_sin)         (using unit normal components)
```

PrusaSlicer uses the formula (Vojtech's triangle area error metric):

```cpp
height = 1.44 * max_surface_deviation * sqrt(n_sin / n_cos)
```

where:
- `n_cos = |normal.z|` (cosine of angle between face normal and Z axis)
- `n_sin = sqrt(normal.x² + normal.y²)` (sine component)
- `max_surface_deviation` is derived from the quality factor (see below)

An additional clamp prevents unrealistically large layer heights on near-vertical surfaces:

```cpp
height = min(max_surface_deviation / 0.184f, computed_height)
```

The constant `0.184` (empirically determined by Wasserfall) represents the volumetric error
coefficient for stacked elliptic extrusion threads.

## Algorithm Walk-through

### Step 1: Mesh Preparation (`SlicingAdaptive::prepare`)

1. The object's triangle mesh is extracted and transformed with the instance matrix.
2. For each triangle face, three values are computed and stored:
   - **z_span**: `(min_z, max_z)` — the Z-range the face occupies.
   - **n_cos**: `|normal.z|` — how "horizontal" the face is (1 = flat top, 0 = vertical wall).
   - **n_sin**: `sqrt(normal.x² + normal.y²)` — how "vertical" the face is.
3. Faces are **sorted by z_span** to allow efficient Z-range queries during layer generation.

### Step 2: Quality Factor → Max Surface Deviation

The user-facing **quality slider** maps to a `quality_factor` in `[0.0, 1.0]`:

| Quality Factor | Meaning | Resulting Layer Height |
|---|---|---|
| 0.0 | Highest quality | `min_layer_height` |
| 0.5 | Standard quality | `layer_height` (default) |
| 1.0 | Lowest quality (fastest) | `max_layer_height` |

The mapping uses linear interpolation (piecewise, through the midpoint):

```cpp
max_surface_deviation = (quality_factor < 0.5)
    ? lerp(delta_min, delta_mid, 2.0 * quality_factor)
    : lerp(delta_max, delta_mid, 2.0 * (1.0 - quality_factor));
```

where `delta_min`, `delta_mid`, `delta_max` correspond to `min_layer_height`,
`layer_height`, and `max_layer_height` respectively.

### Step 3: Layer-by-Layer Height Computation (`SlicingAdaptive::next_layer_height`)

Starting from the top of the previous layer at height `print_z`, the algorithm:

1. **Initializes** `height = max_layer_height`.
2. **Scans all faces** whose Z-span overlaps `print_z` (i.e., the face is intersected by the
   current print plane):
   - For each such face, calls `layer_height_from_slope()` to compute the maximum allowable
     layer height for that face's slope.
   - Takes the **minimum** across all intersecting faces.
3. **Applies min/max clamps**: `height = clamp(height, min_layer_height, max_layer_height)`.
4. **Lookahead correction**: Scans faces whose Z-span *starts* within the proposed layer height.
   If a newly encountered face would require a smaller height, the layer is shrunk accordingly
   (and the process repeats with the lower height). This prevents a layer from being too thick
   and "hiding" a steep feature that begins mid-layer.

### Step 4: Layer Profile Storage

The result is a vector of `(z, height)` pairs — the **layer height profile**:

```
[0.0, h0,  z1, h1,  z2, h2,  ...]
```

This profile is stored in `ModelObject::layer_height_profile` and later consumed by
`generate_object_layers()` to produce the final layer boundaries used during slicing.

## Configuration Parameters

| Parameter | Where Set | Description |
|---|---|---|
| `layer_height` | PrintObjectConfig | Default/target layer height (quality = 0.5) |
| `min_layer_height` | PrintConfig (per nozzle) | Minimum layer height; defaults to 0.07 mm |
| `max_layer_height` | PrintConfig (per nozzle) | Maximum layer height; defaults to 0.75 × nozzle diameter |
| Quality factor | GUI slider | Controls the balance between quality and speed |

## GUI Integration

The adaptive layer height feature is exposed in the **Layer Editing** panel of the 3D view:

1. The user opens the layer height editor and moves the **quality slider**.
2. `GLCanvas3D::adaptive_layer_height_profile(quality_factor)` is called.
3. This calls `layer_height_profile_adaptive(slicing_params, object, quality_factor)`.
4. The resulting profile is stored in `ModelObject::layer_height_profile`.
5. A background re-slice is triggered to apply the new heights.

The user can also **manually edit** the profile by dragging in the layer height editor,
or **smooth** the profile using a Gaussian blur to avoid abrupt height transitions.

## Comparison with Alternative Metrics

Several other layer-height-from-slope formulas exist (shown as comments in `SlicingAdaptive.cpp`):

| Formula | Metric | Notes |
|---|---|---|
| `H = C × tan(α)` | Constant horizontal step distance | Used by Cura; corresponds to topographic contour spacing |
| `H = C × sqrt(tan(α))` | Constant error triangle area (PrusaSlicer) | Good balance of quality and speed |
| `H = C / cos(α)` | Constant perpendicular distance | Old Slic3r metric; used by Perez and Pandey |
| Wasserfall's formula | Volumetric error of elliptic extrusion threads | The original paper's metric |

PrusaSlicer uses the **constant error triangle area** metric with an upper clamp derived from
the empirical roughness constant `0.184`.

## Feature Separation Note

There are two distinct "adaptive" features in PrusaSlicer. Do not confuse them:

1. **Adaptive Layer Height** (this document): Varies *layer thickness* based on model geometry.
   - Implemented in `SlicingAdaptive.hpp/cpp` and `Slicing.cpp`.

2. **Adaptive Cubic Infill** (unrelated): Varies *infill density* using an octree subdivision.
   - Implemented in `Fill/FillAdaptive.hpp/cpp`.
