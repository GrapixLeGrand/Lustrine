#pragma once

#include "Simulation.hpp"
#include "Lustrine_Export.h"

namespace Lustrine {

static inline bool check_index(int i, int min, int max) {
    return (i >= min && i < max);
}

extern void LUSTRINE_EXPORT counting_sort_v2(Simulation* simulation);
extern void find_neighbors_brute_force(Simulation* simulation);
extern void clear_neighbors(Simulation* simulation);
extern void find_neighbors_uniform_grid(Simulation* simulation);
extern void find_neighbors_uniform_grid_v1(Simulation* simulation);
extern void find_neighbors_uniform_grid_v2(Simulation* simulation);
extern void find_neighbors_uniform_grid_v3(Simulation* simulation);
extern void find_neighbors_counting_sort(Simulation* simulation);

};