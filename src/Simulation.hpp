#pragma once



#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
//#include "thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif

//#include "glm/glm.hpp"
#include <vector> 
#include "BulletPhysics.hpp"

namespace Lustrine {

struct Simulation;
typedef float (*W_fun)(const Simulation*, float);
typedef glm::vec3 (*gradW_fun)(const Simulation*, const glm::vec3&);
typedef void (*Simulate_fun)(Simulation*, float);


/**
 * @brief The phyical type of the particles making an object
 */
enum MaterialType {
    SOLID = 0,
    SAND = 1
    
};

/*
A chunk represents particles as a continous array
of positions and thier respective colors.
*/
struct Chunk {
    std::vector<glm::vec3> positions;
    std::vector<glm::vec4> colors;
    glm::vec4 color;
    std::vector<int> cells;

    bool has_one_color_per_particles = false;

    int num_particles;
    //int particlesX;
    //int particlesY;
    //int particlesZ;
    MaterialType type;
};

/*
Represents a grid of particles with each cell
that can contain (or not) a particle. It is
similar to a voxel grid but with particles
*/
struct Grid {

    std::vector<int> cells;
    std::vector<glm::vec4> colors;
    glm::vec4 color;
    bool has_one_color_per_cell;
    int X;
    int Y;
    int Z;
    int num_grid_cells;
    int num_occupied_grid_cells;
    MaterialType type;
    glm::vec3 position;

    //bellow are to be used to parametrise bullet
    bool sparse_solid;//true imply a box per occupied cell, false imply a single box (ignored if gridtype is not solid)
    bool dynamic_solid;//true imply that the collider will be either dynamic or static (ignore if gridtype is not SOLID and if sparse solid is true)

};

/*
Represents the initial parameters of the simulation
The grid array will be given to the simulation to generate
chunks
*/
struct SimulationParameters {

    int X; //domain X
    int Y; //domain Y
    int Z; //domain Z

    float particleRadius;
    float particleDiameter;

};

struct ParticleSink {
    std::vector<bool> state;
    std::vector<float> timers;
    std::vector<float> frequencies;
    std::vector<int> despawned;
    std::vector<std::vector<int>> sink_cells;
    int num_sinks = 0;

    std::vector<int> temp_removal;
};

struct ParticleSource {
    std::vector<Chunk> patterns;
    std::vector<glm::vec3> directions;
    std::vector<float> frequencies; //1/60 means 60 fps spawning
    std::vector<float> timers;
    std::vector<int> capacities;
    std::vector<int> spawned;
    std::vector<bool> source_state;//true=ongoing, false=stop (irrevlantly of the capacity)
    int num_sources = 0;
};

struct WindSystem {
    glm::vec3 direction;
    float magnitude;//read only
    float peak;//max magnitude reached by the wind
    
    float freq_blow;//frequency in seconds of interval between blows
    float timer_blow;//dedicated timer

    float t1;//time to reach peak magnitude
    float t2;//time to reach 0 magniture again
};

struct CountingSortArrays {
    int* counts = nullptr; //size num_grid_cells + 1
    int* particles_unsorted_indices = nullptr;//size num_sand_particles
    int* particles_sorted_indices = nullptr;
};


struct Simulation {

    Simulate_fun simulate_fun = nullptr;
    Bullet::Simulation bullet_physics_simulation;
    W_fun W = nullptr;//pointer to function representing the kernel
    gradW_fun gradW = nullptr;//pointer to function representing the gradient of the kernel

    SimulationParameters parameters_copy;
    int subdivision;
    float particleRadius;
    float particleDiameter;
    float kernelRadius;
    float kernelFactor = 0.5f;

    //stuff for the kernel
    float cubic_kernel_k; //intra parameters to the kernel
    float cubic_kernel_l;

    float domainX = 30.0f; //sizes of the domain
    float domainY = 35.0f;
    float domainZ = 30.0f;

    float rest_density = 24.0; //rest density of fluid
    float mass = 5.0; //mass of each particle

    glm::vec3 gravity = glm::vec3(0, -10.0, 0.0);
    float time_step = 0.01f;

    float relaxation_epsilon = 10.0f;

    //intrisinc parameter to fluid
    float s_corr_dq = 0.5f;
    float s_corr_k = 1.0;
    float s_corr_n = 4;

    float c_xsph = 0.1f;
    float epsilon_vorticity = 0.1f;
    
    glm::vec3* velocities = nullptr;//array containing the velocities of each particles (sand only)
    std::vector<float> lambdas;
    std::vector<std::vector<int>> neighbors; //arrays of neighbor indices
    int gridX, gridY, gridZ; //sizes of the grid
    float cell_size; //size of side length of a single grid cell
    int num_grid_cells; //total amount of grid cells
   
    //uniform grid cells array
    std::vector<std::vector<int>> uniform_gird_cells;
    std::vector<std::vector<int>> uniform_grid_cells_static_saved;//stores the static particles indices (cache them)
    bool computed_static_particles = false;//whether or not we computed the stored indices
    std::vector<std::pair<int, int>> sand_particle_cell_id;
    glm::vec3* velocity_tmp = nullptr;
    glm::vec3* position_star_neighbor_tmp = nullptr;
    glm::vec3* position_neighbor_tmp = nullptr;


    // TODO Everything for the sand goes bellow /////////////////////////////7
    int num_particles = 0; // total amount of particles in the system

    size_t total_allocated = 0;//total (in float) allocated 
    size_t leftover_allocated = 0;//leftover (in float) allocated

    //double list
    glm::vec3* positions = nullptr;//array containing positions of the particles of sand and then space and then solid
    glm::vec3* positions_star = nullptr;
    glm::vec3* positions_solid = nullptr;//array containing positions of the particles (all)

    glm::vec4* colors = nullptr;
    glm::vec4* colors_solid = nullptr;

    int ptr_sand_start = -1;
    int ptr_sand_end = -1;

    int ptr_solid_start = -1;
    int ptr_solid_end = -1;

    int ptr_solid_ordered_start = -1;
    int ptr_solid_ordered_end = -1;

    //We are going to store solid and sand particles not in the same place
    int num_solid_particles;//total amount of particles of solid
    std::vector<Grid> grids_solid;
    std::vector<glm::vec3> grids_initial_positions_solid;
    std::vector<std::pair<int, int>> solid_grid_to_body;//grid index to body index in bullet (WARNING can be -1 bullet index)
    std::vector<std::pair<int, int>> grids_solid_chunk_ptrs;
    std::vector<Chunk> chunks_solid;//all the chunks of solids

    int num_sand_particles;//the total amount of particles of sand
    std::vector<Grid> grids_sand;
    std::vector<glm::vec3> grids_initial_positions_sand;
    std::vector<Chunk> chunks_sand;//all the chunks of sands
    int num_remaining_sand_particles;
    //utils
    glm::vec3* positions_tmp = nullptr;
    bool first_iteration = true;

    //for blowing, attracting sand
    bool attract_flag = false;
    bool blow_flag = false;
    bool* attracted = nullptr;
    bool* attracted_tmp = nullptr;

    float attract_radius = 1.5f;
    float blow_radius = 2.0f;
    float attract_coeff = 1000.0f;
    float blow_coeff = 500.0f;

    ParticleSource* source;
    ParticleSink* sink;
    WindSystem* wind_system;
    CountingSortArrays* counting_sort_arrays;

    float total_time = 0.0f;

};



};