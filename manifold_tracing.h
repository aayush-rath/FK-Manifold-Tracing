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
#include <unordered_map>
#include <queue>

bool edge_intersection(
    Permutahedral_Simplex&s,
    const FK_Triangulation& fk,
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

    int iters = 100;
    double eps = 1e-2 / fk.scale;
    double f0 = sdf(p0);
    double f1 = sdf(p1);

    if (sdf(p0) * sdf(p1) > 1e-4) return false;

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

bool bound_check(
    const FK_Triangulation& fk,
    double* point
) {
    for (int i = 0; i < fk.amb_dim; i++) if (point[i] > 3.14 || point[i] < -3.14) return false;
    return true;
}

void traceManifold(
    const FK_Triangulation& fk,
    double (*sdf)(const double*),
    std::unordered_map<Permutahedral_Simplex, double*, Permutahedral_Simplex_Hash>& Ls,
    double* seed
) {
    std::queue<Permutahedral_Simplex> Q;

    Permutahedral_Simplex initial_simplex = locate_simplex(fk, seed);

    Permutahedral_Simplex initial_edges[MAX_FACES];
    int num_faces = faces(initial_simplex, initial_edges, 1);
    std::cout << num_faces << std::endl;

    for (int i = 0; i < num_faces; i++) {
        Q.push(initial_edges[i]);
        double intersection_point[MAX_D];
        bool h = edge_intersection(initial_edges[i], fk, sdf, intersection_point);
    }

    int intersect_count =  0;

    while (!Q.empty()) {
        Permutahedral_Simplex edge = Q.front();
        Q.pop();

        Permutahedral_Simplex triangles [MAX_COFACES];
        int num_cofaces = cofaces(edge, triangles, 2);
        for (int i = 0; i < num_cofaces; i++) {
            Permutahedral_Simplex new_edges[MAX_FACES];
            int num_edges = faces(triangles[i], new_edges, 1);
            for (int j = 0; j < num_edges; j++) {
                double intersection_point[MAX_D];
                if (!edge_intersection(new_edges[j], fk, sdf, intersection_point)) continue;
                if (Ls.find(new_edges[j]) != Ls.end()) continue;
                if (!bound_check(fk, intersection_point)) continue;
                intersect_count++;
                Q.push(new_edges[j]);
                Ls[new_edges[j]] = intersection_point;
            }
        }

        if (intersect_count > 2000) break;
    }

    std::cout << "Intersection Count: " << intersect_count << std::endl;
}