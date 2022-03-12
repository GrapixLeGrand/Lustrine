#include "Lustrine.hpp"

#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace Lustrine {

constexpr double pi = 3.14159265358979323846;

void init_simulation(const SimulationParameters* parameters, Simulation* simulation, std::vector<Grid> grids, std::vector<glm::vec3> positions) {

    //init_bullet(&simulation->bullet_physics_simulation);

    simulation->domainX = parameters->X;
    simulation->domainY = parameters->Y;
    simulation->domainZ = parameters->Z;

    simulation->chunks.reserve(grids.size());
    for (int i = 0; i < grids.size(); i++) {
        Chunk chunk;
        init_chunk_from_grid(parameters, &chunk, &grids[i], positions[i], grids[i].type);
        simulation->chunks.push_back(chunk);
    }

    std::sort(
        simulation->chunks.begin(),
        simulation->chunks.end(), 
        [](const auto a, const auto b) {
            return a.type < b.type; 
        }
    );

    int num_particles = 0;
    for (int i = 0; i < simulation->chunks.size(); i++) {
        num_particles += simulation->chunks[i].num_particles;
    }

    simulation->num_particles = num_particles;

    simulation->positions = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->positions_star = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->velocities = std::vector<glm::vec3>(simulation->num_particles, {0, 0, 0});
    simulation->lambdas = std::vector<float>(simulation->num_particles, 0);

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

    int current_type = simulation->chunks[0].type;
    simulation->ptr_fluid_start = 0;
    simulation->ptr_static_start = -1;
    simulation->ptr_static_end = simulation->num_particles;
    int offset = 0;

    for (int i = 0; i < simulation->chunks.size(); i++) {
        Chunk& chunk = simulation->chunks[i];
        if (chunk.type != current_type) {
            simulation->ptr_fluid_end = offset;
            simulation->ptr_static_start = offset;
        }
        for (int j = 0; j < simulation->chunks[i].num_particles; j++) {
            simulation->positions[j + offset] = simulation->chunks[i].positions[j];
            simulation->positions_star[j + offset] = simulation->chunks[i].positions[j];
            if (simulation->chunks[i].has_one_color_per_particles) {
                simulation->colors[j + offset] = simulation->chunks[i].colors[j];
            } else {
                simulation->colors[j + offset] = simulation->chunks[i].color;
            }
        }
        offset += simulation->chunks[i].num_particles;
    }

    if (simulation->ptr_static_start == -1) { //no static particles
        simulation->ptr_fluid_end = simulation->num_particles - 1;
        simulation->ptr_static_start = simulation->ptr_static_end;
    }
}


void clean_simulation(Simulation* simulation) {
    clean_bullet(&simulation->bullet_physics_simulation);
}

void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, MaterialType type, glm::vec4 color) {
    
    grid->type = type;
    grid->num_grid_cells = X * Y * Z;
    grid->num_occupied_grid_cells = grid->num_grid_cells;
    grid->X = X;
    grid->Y = Y;
    grid->Z = Z;
    grid->has_one_color_per_cell = false;

    grid->cells = std::vector<bool>(grid->num_grid_cells, true);
    grid->color = color;

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;

}

void init_chunk_from_grid(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, glm::vec3 position, MaterialType type) {
    
    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells;
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, {0, 0, 0});
    if (chunk->has_one_color_per_particles) {
        chunk->colors = std::vector<glm::vec4>(chunk->num_particles, {0, 0, 0, 1});
    } else {
        chunk->color = grid->color;
    }

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;
    glm::vec3 offset = glm::vec3(radius, radius, radius);
    int X = grid->X; int Y = grid->Y; int Z = grid->Z;
    int counter = 0; //because chunk is disorganized we cannot use grid indices

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                
                if (grid->cells[x * Y * Z + y * Z + z] == true) { //a particle present at the current grid cell
                    glm::vec3& particle_position = chunk->positions[counter];

                    particle_position.x = x * diameter;
                    particle_position.y = y * diameter; 
                    particle_position.z = z * diameter; 

                    particle_position += position;

                    if (chunk->has_one_color_per_particles == true) {
                        const glm::vec4& grid_color = grid->colors[x * Y * Z + y * Z + z];
                        glm::vec4& chunk_color = chunk->colors[counter];
                        chunk_color = grid_color;
                    }

                    particle_position += offset;
                    counter++;
                }

            }
        }
    }
}

/*
extern void init_chunk_box(const SimulationParameters* parameters, Chunk* chunk, int X, int Y, int Z, glm::vec3 position, MaterialType type, glm::vec4 chunkColor) {
    
    chunk->type = type;
    chunk->num_particles = X * Y * Z;
    chunk->particlesX = X;
    chunk->particlesY = Y;
    chunk->particlesZ = Z;
    //chunk->has_one_color_per_particles = false;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, {0, 0, 0});
    chunk->colors = std::vector<glm::vec4>(chunk->num_particles, {0, 0, 0, 1});

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;
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

void init_chunk_box(const SimulationParameters* parameters, Chunk* chunk, int X, int Y, int Z) {
    init_chunk_box(parameters, chunk, X, Y, Z, glm::vec3(0.0), MaterialType::FLUID_DYNAMIC, glm::vec4(1.0));
}

/////////

*/


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
    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;
    
    float kernelRadius = simulation->kernelRadius;

    //integration
    for (int i = simulation->ptr_fluid_start; i < simulation->ptr_fluid_end; i++) {
        velocities[i] += simulation->gravity * simulation->mass * dt;
        positions_star[i] = positions[i] + velocities[i] * dt; //prediction
    }

    //find_neighbors_counting_sort(simulation);
    find_neighbors_uniform_grid(simulation);
    //find_neighbors_brute_force(simulation);

    //solve pressure
    for (int i = simulation->ptr_fluid_start; i < simulation->ptr_fluid_end; i++) {

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
    
    for (int i = simulation->ptr_fluid_start; i < simulation->ptr_fluid_end; i++) {

        //equation 13 (applying pressure force correction)
        //pressures_forces[i] = glm::vec3(0.0);
        glm::vec3 pressure_force = glm::vec3(0.0);
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 ij = positions_star[i] - positions_star[neighbors[i][j]];
            pressure_force += (lambdas[i] + lambdas[j] + s_coor(simulation, glm::length(ij))) * simulation->gradW(simulation, ij);
        }

        pressure_force /= simulation->rest_density;
        
        //update prediction
        positions_star[i] += pressure_force;

        //update velocity
        positions_star[i].x = resolve_collision(positions_star[i].x, simulation->particleRadius, X - simulation->particleRadius);
        positions_star[i].y = resolve_collision(positions_star[i].y, simulation->particleRadius, Y - simulation->particleRadius);
        positions_star[i].z = resolve_collision(positions_star[i].z, simulation->particleRadius, Z - simulation->particleRadius);

        velocities[i] = (positions_star[i] - positions[i]) / simulation->time_step;
        positions[i] = positions_star[i];

    }

}

/*
extern void init_chunk_from_grid(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, glm::vec3 position, MaterialType type) {

    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells;
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, {0, 0, 0});
    chunk->colors = std::vector<glm::vec4>(chunk->num_particles, {0, 0, 0, 0});

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;
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

                particle_position.x = x * diameter + position.x;
                particle_position.y = y * diameter + position.y;
                particle_position.z = z * diameter + position.z;

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

}*/
/*
void init_grid_from_magika_voxel(Grid* grid, const std::string& path) {
    init_grid_from_magika_voxel_dont_call_me(grid, path);
}*/

};