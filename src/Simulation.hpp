#pragma once

#include "glm/glm.hpp"
#include <vector> 

namespace Lustrine {

struct Simulation;
typedef float (*W_fun)(const Simulation*, float);
typedef glm::vec3 (*gradW_fun)(const Simulation*, const glm::vec3&);

struct ParticleInCell {
    int cell_id; //cell index of the particle in the grid
    int index; //index of the particle in the unsorted array
};

struct Simulation {
    
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