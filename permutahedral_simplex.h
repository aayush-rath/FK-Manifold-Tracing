#pragma once 

/*
Permutahedral Representation to represent any dimensional simplex
in the Freudenthal-Kuhn Triangulation of the ambient space

Aayush Rath
*/

#include "utils.h"

struct Permutahedral_Simplex {
    int32_t anchor[MAX_D];                                                                              // Minimal Vertex of the simplex representing the location of the simplex
    uint8_t amb_dim;                                                                                    // Ambient space dimension
    uint8_t num_blocks;                                                                                 // Number of ordered partitions in the set. Dimension of the simplex = number of blocks

    uint8_t block_sizes[MAX_D+1];                                                                       // The sizes of each of the partitions
    uint8_t blocks[MAX_D+1][MAX_D+1];                                                                   // The ordered partitions
};

// If two simplex have everything the same then they are equal
bool operator==(
    const Permutahedral_Simplex& a,
    const Permutahedral_Simplex& b
) {
    if (a.amb_dim != b.amb_dim) return false;
    for (int i = 0; i < a.amb_dim; i++) {
        if (a.anchor[i] != b.anchor[i]) return false;
    }

    if (a.num_blocks != b.num_blocks) return false;

    for (int i = 0; i < a.num_blocks; i++) {
        if (a.block_sizes[i] != b.block_sizes[i]) return false;
        for (int j = 0; j < a.block_sizes[i]; j++) {
            if (a.blocks[i][j] != b.blocks[i][j]) return false;
        }
    }

    return true;
}

struct Permutahedral_Simplex_Hash{
    std::size_t operator()(const Permutahedral_Simplex& s) const {
        std::size_t h = 0;

        auto hash_combine = [&](std::size_t v) {
            h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };

        hash_combine(s.amb_dim);
        hash_combine(s.num_blocks);

        for (int i = 0; i < s.amb_dim; i++) {
            hash_combine(std::hash<int32_t>{}(s.anchor[i]));
        }

        for (int i = 0; i < s.num_blocks; i++) {
            hash_combine(s.block_sizes[i]);
            for (int j = 0; j < s.block_sizes[i]; j++) hash_combine(s.blocks[i][j]); 
        }

        return h;
    }
};

// Returns the k-dimensional faces of an l-dimensional simplex
int faces(
    const Permutahedral_Simplex& s,
    Permutahedral_Simplex *k_faces,
    uint8_t k
) {
    const uint8_t l = s.num_blocks - 1;
    if (k > l) return 0;

    uint8_t comb[MAX_D];
    init_combination(comb, k);

    int face_index = 0;

    do {
        Permutahedral_Simplex k_face;
        k_face.amb_dim = s.amb_dim;

        for (int i = 0; i < s.amb_dim; i++) k_face.anchor[i] = s.anchor[i];
        for (int i = 0; i <= k; i++) k_face.block_sizes[i] = 0;

        for (int b = 0; b < comb[0]; b++) {
            for (int j = 0; j < s.block_sizes[b]; j++) {
                int idx = s.blocks[b][j];
                if (idx == s.amb_dim) {
                    for (int d = 0; d < s.amb_dim; ++d)
                        k_face.anchor[d]--;
                } else {
                    k_face.anchor[idx]++;
                }
            }
        }

        k_face.num_blocks = k + 1;

        for (int i = 0; i <= k; i++) {
            if (i < k) {
                k_face.block_sizes[i] = 0;
                for (int j = comb[i]; j < comb[i+1]; j++) {
                    for (int b = 0; b < s.block_sizes[j]; b++) k_face.blocks[i][k_face.block_sizes[i]++] = s.blocks[j][b];
                }
            }

            if (i == k) {
                k_face.block_sizes[k] = 0;
                for (int j = 0; j < comb[0]; j++) {
                    for(int b = 0; b < s.block_sizes[j]; b++) k_face.blocks[k][k_face.block_sizes[k]++] = s.blocks[j][b];
                }
                for (int j = comb[k]; j < l + 1; j++) {
                    for(int b = 0; b < s.block_sizes[j]; b++) k_face.blocks[k][k_face.block_sizes[k]++] = s.blocks[j][b];
                }
            }
        }
        k_faces[face_index++] = k_face;
    } while (next_combination(comb, k, l));

    return face_index;
}