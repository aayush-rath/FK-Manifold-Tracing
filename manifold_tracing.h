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

using Point = std::array<double, MAX_D>;


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

    // std::cout << "Point 1: ";
    // for (int i = 0; i < fk.amb_dim; i++) std::cout << p0[i] << " ";
    // std::cout << std::endl;

    // std::cout << "Point 2: ";
    // for (int i = 0; i < fk.amb_dim; i++) std::cout << p1[i] << " ";
    // std::cout << std::endl;

    int iters = 100;
    double eps = 1e-10 / fk.scale;
    double f0 = sdf(p0);
    double f1 = sdf(p1);

    // std::cout << "SDF Point 1: " << f0 << " SDF Point 2: " << f1 << " Tolerance: " << eps << std::endl;
    //  if (f0 * f1 > 0.0) std::cout << "It doesn't intersect!" << std::endl;

    if (f0 * f1 > 0.0) return false;
    // std::cout << "It intersects!" << std::endl;

    for (int i = 0; i < s.amb_dim; i++) intersection_point[i] = p0[i];

    double f_mid = sdf(p0);
    // std::cout << "Sdf value at mid: " << std::abs(f_mid) << " " << (std::abs(f_mid) > eps? "True" : "False") << std::endl;

    int iter = 0;
    while (std::abs(f_mid) > eps) {
        iter++;

        for (int i = 0; i < s.amb_dim; i++) intersection_point[i] = (f1 * p0[i] - f0 * p1[i])/(f1 - f0);
        f_mid = sdf(intersection_point);
        if (f_mid * f1 > 0) {
            for (int i = 0; i < s.amb_dim; i++) p1[i] = intersection_point[i];
            f1 = sdf(p1);
        } else {
            for (int i = 0; i < s.amb_dim; i++) p0[i] = intersection_point[i];
            f0 = sdf(p0);
        }

        if (iter == iters) break;
    }

    // std::cout << "Calculated Intersection Point: ";
    // for (int i = 0; i < fk.amb_dim; i++) std::cout << intersection_point[i]  << " ";
    // std::cout << std::endl;

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
    std::unordered_map<Permutahedral_Simplex, Point, Permutahedral_Simplex_Hash>& Ls,
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
        Point p;
        for (int i = 0; i < fk.amb_dim; i++) p[i] = intersection_point[i];
        if (h) Ls[initial_edges[i]] = p;
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
                if (Ls.find(new_edges[j]) != Ls.end()) continue;
                if (!edge_intersection(new_edges[j], fk, sdf, intersection_point)) continue;
                // if (!bound_check(fk, intersection_point)) continue;
                intersect_count++;
                Q.push(new_edges[j]);
                // std::cout << "Intersection Point: ";
                // for(int i = 0; i < fk.amb_dim; i++) std::cout << intersection_point[i] << " ";
                // std::cout << std::endl;
                Point p;
                for (int i = 0; i < fk.amb_dim; i++)
                    p[i] = intersection_point[i];

                Ls[new_edges[j]] = p;
                if (intersect_count % 10000 == 0)
                {std::cout << "Intersection Point: ";
                for(int i = 0; i < fk.amb_dim; i++) std::cout << Ls[new_edges[j]][i] << " ";
                std::cout << std::endl;}
            }
        }

        if (intersect_count > 1000000) {
            std::cout << "Size of the queue: " << Q.size() << std::endl;
            break;
        }
    }

    std::cout << "Intersection Count: " << intersect_count << std::endl;
}

void triangulate_surface(
    FK_Triangulation& fk,
    std::unordered_map<Permutahedral_Simplex, Point, Permutahedral_Simplex_Hash>& Ls,
    std::unordered_map<Permutahedral_Simplex, std::vector<Point>, Permutahedral_Simplex_Hash>& Ps
) {
    std::cout << "=== TRIANGULATE SURFACE ===" << std::endl;
    std::cout << "Input edges: " << Ls.size() << std::endl;
    
    int total_assignments = 0;
    
    for (const auto& [edge, point] : Ls) {
        
        // Get all d-dimensional cofaces
        Permutahedral_Simplex d_simplices[MAX_COFACES];
        int num_cofaces = cofaces(edge, d_simplices, fk.amb_dim);
        
        if (total_assignments == 0) {
            std::cout << "First edge has " << num_cofaces << " d-cofaces" << std::endl;
        }

        // std::cout << "Point: ";
        // for (int i = 0; i < fk.amb_dim; i++) std::cout << point[i] << " ";
        // std::cout << std::endl;
        
        // Add this point to all cofaces
        for (int i = 0; i < num_cofaces; i++) {
            Ps[d_simplices[i]].push_back(point);
            total_assignments++;
        }
    }
    
    std::cout << "Total d-simplices: " << Ps.size() << std::endl;
    std::cout << "Total assignments: " << total_assignments << std::endl;
    
    // Check distribution
    int with_d_points = 0;
    for (const auto& [simplex, points] : Ps) {
        if (points.size() == fk.amb_dim) {
            with_d_points++;
        }
    }
    std::cout << "d-simplices with exactly d points: " << with_d_points << std::endl;
    std::cout << "========================" << std::endl;
}