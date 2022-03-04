#pragma once

#include "Simulation.hpp"
#include <string>

namespace Lustrine {
    extern void init_grid_from_magika_voxel(Grid* grid, const std::string& path);
};