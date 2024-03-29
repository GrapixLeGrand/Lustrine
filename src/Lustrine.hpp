#pragma once


#include "Simulation.hpp"
#include "Kernels.hpp"
#include "neighbors/Neighbors.hpp"
#include "Lustrine_Export.h"
#include <vector>
#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
#endif

//

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
    const std::vector<Grid>& grids_solid_arg,
    int subdivision
);


 extern void LUSTRINE_EXPORT init_simulation_extra_parameters(
     const SimulationParameters* parameters,
     Simulation* simulation,
     const std::vector<Grid>& grids_sand_arg,
     const std::vector<Grid>& grids_solid_arg,
     int subdivision,
     float kernel_radius_scale,
     bool with_credits
 );


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
    const std::vector<Grid>& grids_solid_arg
);



/**
 * @brief Clean any allocated and unmanaged memory.
 *
 * @param simulation
 */
extern void LUSTRINE_EXPORT clean_simulation(Simulation* simulation);

/**
 * @brief 
 * 
 * @param chunk 
 * @param grid 
 * @param type 
 * @param cell_size
 * @param subdivision
 */
extern void LUSTRINE_EXPORT init_chunk_from_grid(Chunk* chunk, const Grid* grid, MaterialType type, float cell_size, int subdivision, bool for_bullet);

extern int LUSTRINE_EXPORT add_particle_source(Simulation* simulation, const Grid* pattern, glm::vec3 direction, float freq, int capacity);
extern int LUSTRINE_EXPORT add_particle_sink(Simulation* simulation, glm::vec3 min_pos, glm::vec3 max_pos, float frequency);
extern void LUSTRINE_EXPORT set_source_state(Simulation* simulation, int index, bool state);
extern void LUSTRINE_EXPORT set_sink_state(Simulation* simulation, int index, bool state);

extern int LUSTRINE_EXPORT get_source_spawned(Simulation* simulation, int index);
extern int LUSTRINE_EXPORT get_sink_despawned(Simulation* simulation, int index);

//NYI
extern int LUSTRINE_EXPORT add_particle_sink(Simulation* simulation, const Grid* pattern, float frequency);



/**
 * @brief returns how many particles are in a certain given AABB
 * 
 * @param simulation 
 * @param min 
 * @param max 
 * @param incl_solid 
 * @return int 
 */
extern int LUSTRINE_EXPORT query_cell_num_particles(Simulation* simulation, glm::vec3 min, glm::vec3 max, bool incl_solid);

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
extern void LUSTRINE_EXPORT init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, glm::vec3 position, glm::vec4 color, MaterialType type);
extern void LUSTRINE_EXPORT init_grid_from_magika_voxel(Grid* grid, const std::string& path, glm::vec3 position, MaterialType type);
extern void LUSTRINE_EXPORT init_grid_box_random(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, glm::vec3 position, glm::vec4 color, MaterialType type, float probability);

//extern void fill_grid(Simulation* simulation);
//extern float s_coor(const Simulation* simulation, float rl);

extern void LUSTRINE_EXPORT simulate(Simulation* simulation, float dt); 

/**
 * @brief This function will add a solid chunk to the sim
 * 
 * @param simulation 
 * @return LUSTRINE_EXPORT 
 */
//extern LUSTRINE_EXPORT int add_box(Simulation* simulation, );

//extern LUSTRINE_EXPORT int add_physics_grid_sparse_body();


void update_wind_system(WindSystem* wind_system, float dt);


};