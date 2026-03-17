///|/ Copyright (c) Prusa Research 2016 - 2019 David Kocík @kocikdav, Vojtěch Bubník @bubnikv, Enrico Turri @enricoturri1966
///|/
///|/ PrusaSlicer is released under the terms of the AGPLv3 or higher
///|/
// Based on implementation by @platsch
//
// Adaptive layer height algorithm:
//   For each mesh triangle, the maximum allowable layer height is computed from the
//   triangle's surface slope so that the "stairstepping" surface roughness stays
//   below a user-chosen threshold (the quality factor).  Steeper faces require
//   thinner layers; near-horizontal faces allow thicker layers.
//
// See doc/AdaptiveLayerHeight.md for a full description of the algorithm.

#ifndef slic3r_SlicingAdaptive_hpp_
#define slic3r_SlicingAdaptive_hpp_

#include <stddef.h>
#include <utility>
#include <vector>
#include <cstddef>

#include "Slicing.hpp"
#include "admesh/stl.h"

namespace Slic3r
{

class ModelVolume;
class ModelObject;

class SlicingAdaptive
{
public:
    void  clear();
    void  set_slicing_parameters(SlicingParameters params) { m_slicing_params = params; }
    // Collect all mesh triangles and sort them by their Z-span for fast per-layer queries.
    void  prepare(const ModelObject &object);
    // Return the height of the next layer whose bottom is at print_z.
    // quality_factor is in [0, 1]: 0 = highest quality (thinnest layers),
    //                              0.5 = default layer height,
    //                              1 = lowest quality (thickest layers, fastest print).
    // current_facet is an in/out index into the sorted face list; it is advanced on each
    // call so that successive calls start scanning from where the previous one left off.
	float next_layer_height(const float print_z, float quality_factor, size_t &current_facet);
    // Returns the Z-distance to the next perfectly horizontal facet above z,
    // used to align layer boundaries with flat surface features.
    float horizontal_facet_distance(float z);

    // Per-face data extracted from the mesh and used by the adaptive algorithm.
	struct FaceZ {
        // Minimum and maximum Z coordinates of the triangle's three vertices.
		std::pair<float, float> z_span;
        // |normal.z| – how "horizontal" the face is (1 = flat top, 0 = vertical wall).
		float					n_cos;
        // sqrt(normal.x² + normal.y²) – how "steep" the face is (0 = flat, 1 = vertical).
		float					n_sin;
	};

protected:
	SlicingParameters 		m_slicing_params;

    // Triangles sorted lexicographically by (z_span.first, z_span.second).
	std::vector<FaceZ>		m_faces;
};

}; // namespace Slic3r

#endif /* slic3r_SlicingAdaptive_hpp_ */
