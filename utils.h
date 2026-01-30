#pragma once

/*
Utililty functions to help implement 
the Permutahedral Simplex implementatiton

Aayush Rath
*/


#include <cstdint>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cassert>
#include <cmath>

constexpr int MAX_D = 6;
constexpr int MAX_FACES = 256;
constexpr int MAX_COFACES = 256;

// Get the value of C(n, r) =  number of combinations
int binomial (int n, int r) {
    if (r < 0 || r > n) return 0;
    if (r > n) r = n - r;

    int res = 1;
    for (int i = 0; i < r; i++) res =  res * (n-i) / (i + 1);

    return res;
}

// Initialize the combination with 0, 1, 2, ... k
void init_combination(uint8_t *comb, uint8_t k) {
    for (int i = 0; i <= k ; i++) comb[i] = i;
}

// Get the next combination from the current combination 
bool next_combination(uint8_t *comb, uint8_t k, int l) {
    for (int i = k; i >=0; i--) {
        if (comb[i] < l - (k - i)) {                                                                    // Check if the max possible value is reached
            comb[i]++;                                                                                  // Increment by 1

            for (int j = i+1; j <= k; j++) comb[j] = comb[j-1] + 1;                                     // Reset all the elements right to the current element to the lowest possible value
            return true;
        }
    }
    return false;
}

// Initialize the interger composition
bool walsh_init(uint8_t k, uint8_t l, const uint8_t *bounds, uint8_t* a) {
    uint8_t remaining = l - k;

    for (uint8_t i = 0; i <= k; i++) {
        uint8_t v = remaining;
        if (v > bounds[i]-1) v = bounds[i] - 1;                                                         // Set the values to their bounds
        a[i] = v;                                                                                       // If no more remains the set the value to be the remaining value
        remaining -= v;
    }

    return (remaining == 0);
}

// Give the next integer composition
bool walsh_next(uint8_t k, const uint8_t *bounds, uint8_t *a) {
    for (int i = 0; i < k; i++) {
        if (a[i] == 0) continue;                                                                        // Find the leftmost value that can give something

        for (int j = i + 1; j <= k; j++) {                                                              // Find the next index on right that can give something
            if (a[j] < bounds[j] - 1) {
                a[i]--;
                a[j]++;                                                                                 // Take 1 from a[i] and give to a[j]

                int s = 0;
                for (int t = 0; t < j; t++) {
                    s += a[t];
                    a[t] = 0;
                }                                                                                       // Collect all of the stuff from the left side including the current i-th element

                for (int t = 0; t < j && s > 0; t++) {                                                  // Redistribute the collected stuff including the current element i-th element
                    int v = s;
                    if (v > bounds[t] - 1) v = bounds[t] - 1;
                    a[t] = (uint8_t)v;
                    s -= v;
                }

                return true;
            }
        }
    }

    return false;
}

// Restricted Growth Strings for creating partitions
bool rgs_init(uint8_t omega_i, uint8_t a_i, uint8_t *rgs, uint8_t *max_val) {
    if (omega_i == 0 || a_i == 0 || a_i > omega_i) return false;

    for (uint8_t i = 0; i <= omega_i - a_i; i++) rgs[i] = 0;                                                    // Set all the intial l-k elements into the 0th partition
    for (uint8_t i = omega_i - a_i + 1, j = 1; i < omega_i; i++, j++) rgs[i] = j;                                    // Set all the left out one to the next increment partition since 0

    max_val[0] = 255;
    for (uint8_t i = 1; i <= omega_i + 1; i++) {
        max_val[i] = rgs[i - 1] + 1;
    }

    return true;
}

// Get the next RGS partition
bool rgs_next(uint8_t omega_i, uint8_t a_i, uint8_t *rgs, uint8_t *max_val) {
    if (a_i <= 1) return false;

    int i = omega_i - 1;
    while (i > 0 && (rgs[i] + 1 > max_val[i] || rgs[i] + 1 >= a_i)) i--;                                  // Find the rightmost element in the rgs that can be incremented

    if (i == 0) return false;                                                                           // No more partitions left
    
    rgs[i]++;
    uint8_t mm = max_val[i];
    mm += (rgs[i] >= mm);
    max_val[i+1] = mm;

    while (++i < omega_i + 1) {
        rgs[i] = 0;
        max_val[i + 1] = mm;
    }

    uint8_t p = a_i + 1;
    if (mm < p) {
        do {
            max_val[i] = p;
            i--;
            p--;
            rgs[i] = p;
        } while (max_val[i] < p);
    }
    
    return true;
}