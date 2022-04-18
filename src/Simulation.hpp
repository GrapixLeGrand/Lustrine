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

    std::vector<bool> cells;
    std::vector<glm::vec4> colors;
    glm::vec4 color;
    bool has_one_color_per_cell;
    int X;
    int Y;
    int Z;
    int num_grid_cells;
    int num_occupied_grid_cells;
    MaterialType type;

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

    float particleRadius = 0.5f;
    float particleDiameter = 1.0f;

};

struct Simulation {

    Bullet::Simulation bullet_physics_simulation;
    W_fun W = nullptr;//pointer to function representing the kernel
    gradW_fun gradW = nullptr;//pointer to function representing the gradient of the kernel

    float particleRadius = 0.5;
    float particleDiameter = 2 * particleRadius;
    float kernelRadius = 3.1f * particleRadius;
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
    
    std::vector<glm::vec3> velocities;//array containing the velocities of each particles (sand only)
    std::vector<float> lambdas;
    std::vector<std::vector<int>> neighbors; //arrays of neighbor indices
    int gridX, gridY, gridZ; //sizes of the grid
    float cell_size; //size of side length of a single grid cell
    int num_grid_cells; //total amount of grid cells
    
    //convignent arrays for the grid based methods
    std::vector<std::pair<int, int>> particle_cell_index_to_index;
    std::vector<std::pair<int, int>> cell_indices;
    std::vector<glm::vec3> positions_star_copy;

    //uniform grid cells array
    std::vector<std::vector<int>> uniform_gird_cells;

    //counting sort arrays
    std::vector<int> counts;
    std::vector<int> counting_sort_sorted_indices;

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

    //utils
    glm::vec3* positions_tmp = nullptr;

};



};