#pragma once

#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
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
 * 
 */
enum MaterialType {
    FLUID_DYNAMIC = 0,
    SOLID_STATIC = 1,
    SAND_DYNAMIC = 2
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
    int particlesX;
    int particlesY;
    int particlesZ;
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

    BulletPhyicsSimulation bullet_physics_simulation;

    std::vector<Chunk> chunks;

    int ptr_fluid_start = 0; //start positions of fluid particles
    int ptr_fluid_end = 0;  //end positions of fluid particles

    int ptr_static_start = 0; //start position of solid particles
    int ptr_static_end = 0;   //end position of soldid particles
    
    W_fun W = nullptr; //pointer to function representing the kernel
    gradW_fun gradW = nullptr; //pointer to function representing the gradient of the kernel

    float particleRadius = 0.5;
    float particleDiameter = 2 * particleRadius;
    float kernelRadius = 3.1f * particleRadius;
    float kernelFactor = 0.5f;

    //stuff for the kernel
    float cubic_kernel_k; //intra parameters to the kernel
    float cubic_kernel_l;

    int num_particles; //the total amount of particles

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
    
    std::vector<glm::vec3> positions; //array containing positions of the particles
    std::vector<glm::vec3> positions_star; //array containing the prediction of the positions of the particles
    std::vector<glm::vec3> velocities; //array containing the velocities of each particles

    std::vector<float> lambdas;
    //std::vector<glm::vec3> vorticities; unused for now
    std::vector<std::vector<int>> neighbors; //arrays of neighbor indices
    std::vector<glm::vec4> colors; //per particle colors

    
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

    // TODO Everything for the sand goes bellow
    
};



};