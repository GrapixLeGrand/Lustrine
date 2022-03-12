#pragma once

#include "Simulation.hpp"
#include "Lustrine_Export.h"
#include <string>

namespace Lustrine {
    extern void init_grid_from_magika_voxel_dont_call_me(Grid* grid, const std::string& path, MaterialType type);
    extern LUSTRINE_EXPORT std::string read_vox_scene_json(const uint8_t *buffer, int64_t size);
};