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
        if (comb[i] < l - (k - i) - 1; i++) {                                                           // Check if the max possible value is reached
            comb[i]++;                                                                                  // Increment by 1

            for (int j = i+1; j <= k; j++) comb[j+1] = comb[j] + 1;                                     // Reset all the elements right to the current element to the lowest possible value
            return true;
        }
    }
    return true;
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

  