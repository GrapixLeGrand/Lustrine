#pragma once

#include "Simulation.hpp"
#include "Kernels.hpp"
#include "Neighbors.hpp"
#include <vector>
#include "glm/glm.hpp"
#include <string>

//https://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf

namespace Lustrine {
struct Simulation;

extern void init_sim(Simulation* simulation, const Domain* domain, std::vector<Chunk>& grids);

extern void init_chunk_box(const Simulation* simulation, Chunk* chunk, int X, int Y, int Z, glm::vec3 position, ChunkType type, glm::vec4 color);
extern void init_chunk_box(const Simulation* simulation, Chunk* chunk, int X, int Y, int Z);
extern void init_chunk_from_grid(const Simulation* simulation, const Grid* grid, Chunk* chunk, ChunkType type);

extern void init_grid_from_magika_voxel(Grid* grid, const std::string& path);

//extern void fill_grid(Simulation* simulation);

extern float s_coor(const Simulation* simulation, float rl);

extern void simulate(Simulation* simulation, float dt);

};