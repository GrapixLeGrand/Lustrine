#pragma once


#include "Simulation.hpp"
#include "Kernels.hpp"
#include "Neighbors.hpp"
#include "Lustrine_Export.h"
#include <vector>
#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
#endif

#include <string>

//https://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf

namespace Lustrine {
struct Simulation;

/**
 * @brief Initialize the simulation from the parameters in Simulation parameters.
 * 
 * @param parameters 
 * @param simulation 
 * @param grids, the allocated grids
 * @param the positions of the grids
 */
//extern void LUSTRINE_EXPORT init_simulation(const SimulationParameters* parameters, Simulation* simulation, std::vector<Grid> grids, std::vector<glm::vec3> positions);

/**
 * @brief init simulation with player box shape
 * 
 * @param parameters 
 * @param simulation 
 * @param grids 
 * @param positions 
 * @param player_grid 
 * @param player_position 
 */
extern void LUSTRINE_EXPORT init_simulation(
    const SimulationParameters* parameters, 
    Simulation* simulation, 
    const std::vector<Grid>& grids_sand_arg, 
    const std::vector<glm::vec3>& grids_sand_positions_arg, 
    const std::vector<Grid>& grids_solid_arg, 
    const std::vector<glm::vec3>& grids_solid_positions_arg
);

/**
 * @brief Clean any allocated and unmanaged memory.
 *
 * @param simulation
 */
extern void clean_simulation(Simulation* simulation);

/**
 * @brief 
 * 
 * @param parameters 
 * @param chunk 
 * @param grid 
 * @param position 
 * @param type 
 */
extern void init_chunk_from_grid(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, glm::vec3 position, MaterialType type);

/**
 * @brief Init a grid full of partlices. A box full of particles with X * Y * Z particles.
 * 
 * @param parameters 
 * @param grid 
 * @param X 
 * @param Y 
 * @param Z 
 * @param type 
 * @param color 
 */
extern LUSTRINE_EXPORT void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, MaterialType type, glm::vec4 color);

extern LUSTRINE_EXPORT void init_grid_from_magika_voxel(Grid* grid, const std::string& path, MaterialType type);

//extern void fill_grid(Simulation* simulation);
//extern float s_coor(const Simulation* simulation, float rl);

extern LUSTRINE_EXPORT void simulate(Simulation* simulation, float dt);


//extern LUSTRINE_EXPORT int add_physics_grid_sparse_body();



};