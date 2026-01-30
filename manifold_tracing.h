#pragma once

/*
Manifold Tracing implementation for d-1 manifolds using the Permutahedral 
Representation and finding intersection with 1-simplex

Aayush Rath
*/

#include "utils.h"
#include "permutahedral_simplex.h"
#include "fk_triangulation.h"
#include <cassert>

bool edge_intersection(
    Permutahedral_Simplex&s,
    FK_Triangulation& fk,
    double (*sdf)(const double*),
    double* intersection_point
) {
    assert(s.num_blocks == 2);

    // Endpoints of the edge as FK coordinates
    int32_t v0[MAX_D];
    for (int i = 0; i < s.amb_dim; ++i) v0[i] = s.anchor[i];
    
    int32_t v1[MAX_D];
    for (int i = 0; i < s.amb_dim; ++i) v1[i] = s.anchor[i];

    for (int j = 0; j < s.block_sizes[0]; ++j) {
        int idx = s.blocks[0][j];
        if (idx == s.amb_dim) {
            for (int k = 0; k < s.amb_dim; ++k) v1[k]++;
        } else {
            v1[idx]++;
        }
    }

    // Cartesian coordinates
    double p0[MAX_D], p1[MAX_D];
    fk.cartesian_coordinates(v0, p0);
    fk.cartesian_coordinates(v1, p1);

    if (sdf(p0) * sdf(p1) > 0.0) return false;

    int iters = 100;
    double eps = 1e-2 / fk.scale;
    double f0 = sdf(p0);
    double f1 = sdf(p1);

    int iter = 0;
    while (abs(f0 - f1) > eps) {
        iter++;

        for (int i = 0; i < s.amb_dim; i++) intersection_point[i] = (f1 * p0[i] - f0 * p1[i])/(f1 - f0);
        if (sdf(intersection_point) * f0 > 0) {
            for (int i = 0; i < s.amb_dim; i++) p1[i] = intersection_point[i];
            f1 = sdf(p1);
        } else {
            for (int i = 0; i < s.amb_dim; i++) p0[i] = intersection_point[i];
            f0 = sdf(p0);
        }

        if (iter == iters) break;
    }

    return true;
}