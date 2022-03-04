#include "Lustrine.hpp"

#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace Lustrine {


constexpr double pi = 3.14159265358979323846;

void init_sim(Simulation* simulation, const Domain* domain, std::vector<Chunk>& chunks) {

    std::sort(
        chunks.begin(),
        chunks.end(), 
        [](const auto a, const auto b) {
            return a.type < b.type; 
        }
    );

    int num_particles = 0;

    for (int i = 0; i < chunks.size(); i++) {
        num_particles += chunks[i].num_particles;
    }

    simulation->num_particles = num_particles; //particlesX * particlesY * particlesZ;

    //simulation->max_neighbors = simulation->num_particles - 1;

    simulation->positions = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->positions_star = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->velocities = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->pressures_forces = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});

    simulation->densities = std::vector<float>(simulation->num_particles, 0);
    simulation->lambdas = std::vector<float>(simulation->num_particles, 0);

    simulation->vorticities = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});

    simulation->neighbors = std::vector<std::vector<int>>(simulation->num_particles);
    simulation->colors = std::vector<glm::vec4>(simulation->num_particles, {0, 0, 0, 0});

    simulation->W = cubic_kernel;
    simulation->gradW = cubic_kernel_grad;

    //for the kernel
    float h3 = std::pow(simulation->kernelRadius, 3);
    simulation->cubic_kernel_k = 8.0 / (pi * h3);
    simulation->cubic_kernel_l = 48.0 / (pi * h3);

    //sorting neighbor strategy with grid
    simulation->cell_size = 1.0f * simulation->kernelRadius;
    simulation->gridX = (int) (simulation->domainX / simulation->cell_size) + 1; //plus 1 because its a full range
    simulation->gridY = (int) (simulation->domainY / simulation->cell_size) + 1;
    simulation->gridZ = (int) (simulation->domainZ / simulation->cell_size) + 1;

    simulation->num_grid_cells = (simulation->gridX) * (simulation->gridY) * (simulation->gridZ);
    simulation->particle_cell_index_to_index = std::vector<std::pair<int, int>>(simulation->num_particles, std::make_pair(0, 0));
    simulation->positions_star_copy = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->cell_indices = std::vector<std::pair<int, int>>(simulation->num_grid_cells, std::make_pair(0, 0));

    //uniform grid
    simulation->uniform_gird_cells = std::vector<std::vector<int>> (simulation->num_grid_cells, std::vector<int>{});

    //counting sort
    simulation->counts = std::vector<int>(simulation->num_grid_cells + 1, 0);
    simulation->counting_sort_sorted_indices = std::vector<int>(simulation->num_particles, 0);

    simulation->times = std::vector<double>(10, 0);

    int current_type = chunks[0].type;
    simulation->ptr_fluid_dynamic_start = 0;
    simulation->ptr_fluid_static_start = -1;
    simulation->ptr_fluid_static_end = simulation->num_particles;
    int offset = 0;

    for (int i = 0; i < chunks.size(); i++) {
        Chunk& chunk = chunks[i];
        if (chunk.type != current_type) {
            simulation->ptr_fluid_dynamic_end = offset;
            simulation->ptr_fluid_static_start = offset;
        }
        for (int j = 0; j < chunks[i].num_particles; j++) {
            simulation->positions[j + offset] = chunks[i].positions[j];
            simulation->positions_star[j + offset] = chunks[i].positions[j];
            simulation->colors[j + offset] = chunks[i].colors[j];
        }
        offset += chunks[i].num_particles;
    }

    if (simulation->ptr_fluid_static_start == -1) { //no static particles
        simulation->ptr_fluid_dynamic_end = simulation->num_particles - 1;
        simulation->ptr_fluid_static_start = simulation->ptr_fluid_static_end;
    }

}

extern void init_chunk_box(const Simulation* simulation, Chunk* chunk, int X, int Y, int Z, glm::vec3 position, ChunkType type, glm::vec4 chunkColor) {
    
    chunk->type = type;
    chunk->num_particles = X * Y * Z;
    chunk->particlesX = X;
    chunk->particlesY = Y;
    chunk->particlesZ = Z;
    //chunk->has_one_color_per_particles = false;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, {0, 0, 0});
    chunk->colors = std::vector<glm::vec4>(chunk->num_particles, {0, 0, 0, 1});

    const float diameter = simulation->particleDiameter;
    const float radius = simulation->particleRadius;
    glm::vec3 offset = glm::vec3(radius, radius, radius);

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                
                glm::vec3& particle_position = chunk->positions[x * Y * Z + y * Z + z];
                glm::vec4& color = chunk->colors[x * Y * Z + y * Z + z];

                particle_position.x = x * diameter;
                particle_position.y = y * diameter; 
                particle_position.z = z * diameter; 

                particle_position += position;

                if (chunk->has_one_color_per_particles == true) {
                    color.r = particle_position.x / X;
                    color.g = particle_position.y / Y;
                    color.b = particle_position.z / Z;
                    color.a = 1.0f;
                } else {
                    color = chunkColor;
                }

                particle_position += offset;

            }
        }
    }
}

void init_chunk_box(const Simulation* simulation, Chunk* chunk, int X, int Y, int Z) {
    init_chunk_box(simulation, chunk, X, Y, Z, glm::vec3(0.0), ChunkType::FLUID_DYNAMIC, glm::vec4(1.0));
}

/////////



float s_coor(const Simulation* simulation, float rl) {
    return - simulation->s_corr_k * std::pow(simulation->W(simulation, rl) / simulation->W(simulation, simulation->s_corr_dq), simulation->s_corr_n);
}

float epsilon_collision = 0.01;

float resolve_collision(float value, float min, float max) {
    
    if (value <= min) {
        return epsilon_collision;
    }

    if (value > max) {
        return max - epsilon_collision;
    }

    return value;
}

void simulate(Simulation* simulation, float dt) {

    dt = glm::clamp(dt, 0.001f, 0.01f);
    simulation->time_step = dt;

    int n = simulation->num_particles;
    int X = simulation->domainX;
    int Y = simulation->domainY;
    int Z = simulation->domainZ;

    std::vector<glm::vec3>& positions = simulation->positions;
    std::vector<glm::vec3>& positions_star = simulation->positions_star;
    std::vector<glm::vec3>& velocities = simulation->velocities;
    std::vector<glm::vec3>& pressures_forces = simulation->pressures_forces;

    std::vector<float>& densities = simulation->densities;
    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;
    
    float kernelRadius = simulation->kernelRadius;
    //float dt = simulation->time_step;

    //integration
    for (int i = simulation->ptr_fluid_dynamic_start; i < simulation->ptr_fluid_dynamic_end; i++) {
        velocities[i] += simulation->gravity * simulation->mass * dt;
        positions_star[i] = positions[i] + velocities[i] * dt; //prediction
    }

    //find_neighbors_counting_sort(simulation);
    find_neighbors_uniform_grid(simulation);
    //find_neighbors_brute_force(simulation);

    //solve pressure
    for (int i = simulation->ptr_fluid_dynamic_start; i < simulation->ptr_fluid_dynamic_end; i++) {

        float densitiy = 0.0;
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 ij = positions_star[i] - positions_star[neighbors[i][j]];
            float len = glm::length(ij);
            densitiy += simulation->mass * simulation->W(simulation, len);
        }
        densitiy += simulation->mass * simulation->W(simulation, 0.0);

        //equation 1
        float constraint_i = (densitiy / simulation->rest_density) - 1.0;
        float constraint_gradient_sum = 0.0;
        glm::vec3 grad_current_p = glm::vec3(0.0);

        //equation 8
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 temp = positions_star[i] - positions_star[neighbors[i][j]];
            glm::vec3 neighbor_grad = - (simulation->mass / simulation->rest_density) * simulation->gradW(simulation, temp);
            constraint_gradient_sum += glm::dot(neighbor_grad, neighbor_grad);
            grad_current_p -= neighbor_grad;
        }

        constraint_gradient_sum += glm::dot(grad_current_p, grad_current_p);

        lambdas[i] = 0.0;
        if (constraint_gradient_sum > 0.0) {
            lambdas[i] = - constraint_i / (constraint_gradient_sum + simulation->relaxation_epsilon);
        }

    }
    
    for (int i = simulation->ptr_fluid_dynamic_start; i < simulation->ptr_fluid_dynamic_end; i++) {

        //equation 13 (applying pressure force correction)
        pressures_forces[i] = glm::vec3(0.0);
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 ij = positions_star[i] - positions_star[neighbors[i][j]];
            pressures_forces[i] += (lambdas[i] + lambdas[j] + s_coor(simulation, glm::length(ij))) * simulation->gradW(simulation, ij);
        }

        pressures_forces[i] /= simulation->rest_density;
        
        //update prediction
        positions_star[i] += pressures_forces[i];

        //update velocity
        positions_star[i].x = resolve_collision(positions_star[i].x, simulation->particleRadius, X - simulation->particleRadius);
        positions_star[i].y = resolve_collision(positions_star[i].y, simulation->particleRadius, Y - simulation->particleRadius);
        positions_star[i].z = resolve_collision(positions_star[i].z, simulation->particleRadius, Z - simulation->particleRadius);

        velocities[i] = (positions_star[i] - positions[i]) / simulation->time_step;
        positions[i] = positions_star[i];

    }

}

extern void init_chunk_from_grid(const Simulation* simulation, const Grid* grid, Chunk* chunk, ChunkType type) {

    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells;
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, {0, 0, 0});
    chunk->colors = std::vector<glm::vec4>(chunk->num_particles, {0, 0, 0, 1});

    const float diameter = simulation->particleDiameter;
    const float radius = simulation->particleRadius;
    glm::vec3 offset = glm::vec3(radius, radius, radius);

    int counter = 0;
    for (int x = 0; x < grid->X; x++) {
        for (int y = 0; y < grid->Y; y++) {
            for (int z = 0; z < grid->Z; z++) {
                
                int index = x * grid->Y * grid->Z + y * grid->Z + z;

                if (grid->cells[index] == false) { //no voxel here
                    continue;
                }

                glm::vec3& particle_position = chunk->positions[counter];

                particle_position.x = x * diameter;
                particle_position.y = y * diameter;
                particle_position.z = z * diameter;

                const glm::vec4& voxel_color = grid->colors[index];
                glm::vec4& chunk_color = chunk->colors[counter];

                //particle_position += position;

                if (chunk->has_one_color_per_particles == true) {
                    chunk_color.r = voxel_color.r;
                    chunk_color.g = voxel_color.g;
                    chunk_color.b = voxel_color.b;
                    chunk_color.a = voxel_color.a;
                } else {
                    //color = chunkColor;
                }

                
                counter++;
            } //end for z
        } //end for y
    } //end for x

}



};