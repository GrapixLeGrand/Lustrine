#pragma once

#include "glm/glm.hpp"
#include <vector> 

namespace Lustrine {

struct Simulation;
typedef float (*W_fun)(const Simulation*, float);
typedef glm::vec3 (*gradW_fun)(const Simulation*, const glm::vec3&);

enum ChunkType {
    FLUID_DYNAMIC = 0,
    FLUID_STATIC = 1
};

struct Domain {
    float X = 0;
    float Y = 0;
    float Z = 0;
};


struct Chunk {
    std::vector<glm::vec3> positions;
    std::vector<glm::vec4> colors;
    glm::vec3 color;

    bool has_one_color_per_particles = false;

    int num_particles;
    int particlesX;
    int particlesY;
    int particlesZ;
    ChunkType type;
};

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
};

struct Simulation {

    int ptr_fluid_dynamic_start = 0;
    int ptr_fluid_dynamic_end = 0;

    int ptr_fluid_static_start = 0;
    int ptr_fluid_static_end = 0;
    
    W_fun W = nullptr;
    gradW_fun gradW = nullptr;

    float particleRadius = 0.5;
    float particleDiameter = 2 * particleRadius;
    float kernelRadius = 3.1f * particleRadius;
    float kernelFactor = 0.5f;

    //stuff for the kernel
    float cubic_kernel_k;
    float cubic_kernel_l;

    int max_neighbors;

    int particlesX;
    int particlesY;
    int particlesZ;
    int num_particles;

    float domainX = 30.0f;
    float domainY = 35.0f;
    float domainZ = 30.0f;

    float rest_density = 24.0;
    float mass = 5.0;

    glm::vec3 gravity = glm::vec3(0, -10.0, 0.0);
    float time_step = 0.01;
    float steps = 4;

    float relaxation_epsilon = 10.0f;

    float s_corr_dq = 0.5f;
    float s_corr_k = 1.0;
    float s_corr_n = 4;

    float c_xsph = 0.1;
    float epsilon_vorticity = 0.1;
    
    std::vector<glm::vec3> positions;
    std::vector<glm::vec3> positions_star;
    std::vector<glm::vec3> velocities;
    std::vector<glm::vec3> pressures_forces;

    std::vector<float> densities;
    std::vector<float> lambdas;
    std::vector<glm::vec3> vorticities;
    std::vector<std::vector<int>> neighbors;
    std::vector<glm::vec4> colors;

    
    int gridX, gridY, gridZ;
    float cell_size;
    int num_grid_cells;
    
    std::vector<std::pair<int, int>> particle_cell_index_to_index;
    std::vector<std::pair<int, int>> cell_indices;

    std::vector<glm::vec3> positions_star_copy;

    //uniform grid cells array
    std::vector<std::vector<int>> uniform_gird_cells;

    //counting sort arrays
    std::vector<int> counts;
    std::vector<int> counting_sort_sorted_indices;

    std::vector<double> times;
};



};