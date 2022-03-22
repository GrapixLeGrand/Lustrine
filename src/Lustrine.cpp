#include "Lustrine.hpp"

#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace Lustrine {

constexpr double pi = 3.14159265358979323846;


void add_physics_static_grid(Simulation* simulation, const Grid* grid) {

    for (int x = 0; x < grid->X; x++) {
        for (int y = 0; y < grid->Y; y++) {
            for (int z = 0; z < grid->Z; z++) {
                int index = grid->Y * grid->Z * x + grid->Z * y + z;
                if (grid->cells[index] == true) {
                    glm::vec3 position = {x, y, z};
                    if (grid->has_one_color_per_cell) {
                       add_box(&simulation->bullet_physics_simulation, position, false);     
                    } else {
                        add_box(&simulation->bullet_physics_simulation, position, false);
                    }
                }
            }
        }
    }

}


void init_simulation(const SimulationParameters* parameters, Simulation* simulation, std::vector<Grid> grids, std::vector<glm::vec3> positions) {
    //init_simulation(parameters, simulation, grids, positions, nullptr, {0, 0, 0});
    assert(false);
}

void init_simulation(
    const SimulationParameters* parameters, 
    Simulation* simulation, 
    const std::vector<Grid>& grids_sand_arg, 
    const std::vector<glm::vec3>& grids_sand_positions_arg, 
    const std::vector<Grid>& grids_solid_arg, 
    const std::vector<glm::vec3>& grids_solid_positions_arg) {

    init_bullet(&simulation->bullet_physics_simulation);

    simulation->domainX = parameters->X;
    simulation->domainY = parameters->Y;
    simulation->domainZ = parameters->Z;

    int first_guess_allocation = simulation->domainX * simulation->domainY * simulation->domainZ;
    assert(first_guess_allocation > 0);
    int initial_allocated_particles_num = 0;

    for (int i = 0; i < grids_sand_arg.size(); i++) {
        initial_allocated_particles_num += grids_sand_arg[i].num_occupied_grid_cells;
    }

    simulation->num_sand_particles = initial_allocated_particles_num;

    for (int i = 0; i < grids_solid_arg.size(); i++) {
        initial_allocated_particles_num += grids_solid_arg[i].num_occupied_grid_cells;
    }

    simulation->num_solid_particles = initial_allocated_particles_num - simulation->num_sand_particles;

    simulation->total_allocated = (size_t) std::max(initial_allocated_particles_num, first_guess_allocation);
    simulation->leftover_allocated = simulation->total_allocated - initial_allocated_particles_num;

    simulation->positions = new glm::vec3[simulation->total_allocated];
    simulation->positions_star = new glm::vec3[simulation->total_allocated];
    simulation->colors = new glm::vec4[simulation->total_allocated];

    memset(simulation->positions, 0, simulation->total_allocated * 4 * 3);
    memset(simulation->positions_star, 0, simulation->total_allocated * 4 * 3);
    memset(simulation->colors, 0, simulation->total_allocated * 4 * 4);


    { //sand
        
        simulation->ptr_sand_start = 0;
        simulation->ptr_sand_end = 0;

        simulation->grids_sand = grids_sand_arg;
        simulation->grids_initial_positions_sand = grids_sand_positions_arg;
        simulation->chunks_sand.reserve(simulation->grids_sand.size());

        for (int i = 0; i < simulation->grids_sand.size(); i++) {
            Chunk chunk;
            assert(simulation->grids_sand[i].type == SAND);
            init_chunk_from_grid(parameters, &chunk, &simulation->grids_sand[i], simulation->grids_initial_positions_sand[i], SAND);
            simulation->chunks_sand.push_back(chunk);

            for (int j = 0; j < chunk.num_particles; j++) {
                simulation->positions[simulation->ptr_sand_end] = chunk.positions[j];
                simulation->positions_star[simulation->ptr_sand_end] = chunk.positions[j];
                if (chunk.has_one_color_per_particles == true) {
                    simulation->colors[simulation->ptr_sand_end] = chunk.colors[j];
                } else {
                    simulation->colors[simulation->ptr_sand_end] = chunk.color;
                }
                simulation->ptr_sand_end++;
            }

        }
    } // end sand

    {  //solids
        
        simulation->ptr_solid_start = simulation->total_allocated - simulation->num_solid_particles; //simulation->total_allocated - 1; //head is going down so list is in order
        simulation->ptr_solid_end = simulation->ptr_solid_start;

        simulation->grids_solid = grids_solid_arg;
        simulation->grids_initial_positions_solid = grids_solid_positions_arg;
        simulation->chunks_solid.reserve(simulation->grids_solid.size());

        for (int i = 0; i < simulation->grids_solid.size(); i++) {
            Chunk chunk;
            assert(simulation->grids_solid[i].type == SOLID);
            init_chunk_from_grid(parameters, &chunk, &simulation->grids_solid[i], simulation->grids_initial_positions_solid[i], SOLID);
            simulation->chunks_solid.push_back(chunk);

            for (int j = 0; j < chunk.num_particles; j++) {
                simulation->positions[simulation->ptr_solid_end] = chunk.positions[j];
                simulation->positions_star[simulation->ptr_solid_end] = chunk.positions[j];
                if (chunk.has_one_color_per_particles == true) {
                    simulation->colors[simulation->ptr_solid_end] = chunk.colors[j];
                } else {
                    simulation->colors[simulation->ptr_solid_end] = chunk.color;
                }
                simulation->ptr_solid_end++;
            }

        }

        assert(simulation->ptr_solid_end == simulation->total_allocated);

    } // end solids

    simulation->velocities = std::vector<glm::vec3>(simulation->num_sand_particles, {0.0f, 0.0f, 0.0f});
    simulation->lambdas = std::vector<float>(simulation->num_sand_particles, 0.0f);
    simulation->neighbors = std::vector<std::vector<int>>(simulation->num_sand_particles, std::vector<int>{});

    simulation->W = cubic_kernel;
    simulation->gradW = cubic_kernel_grad;

    //for the kernel
    float h3 = std::pow(simulation->kernelRadius, 3);
    simulation->cubic_kernel_k = 8.0f / (pi * h3);
    simulation->cubic_kernel_l = 48.0f / (pi * h3);

    //sorting neighbor strategy with grid
    simulation->cell_size = 1.0f * simulation->kernelRadius;
    simulation->gridX = (int) (simulation->domainX / simulation->cell_size) + 1; //plus 1 because its a full range
    simulation->gridY = (int) (simulation->domainY / simulation->cell_size) + 1;
    simulation->gridZ = (int) (simulation->domainZ / simulation->cell_size) + 1;

    simulation->num_grid_cells = (simulation->gridX) * (simulation->gridY) * (simulation->gridZ);
    simulation->particle_cell_index_to_index = std::vector<std::pair<int, int>>(simulation->num_sand_particles, std::make_pair(0, 0));
    simulation->positions_star_copy = std::vector<glm::vec3>(simulation->num_sand_particles, {0, 0, 0});
    simulation->cell_indices = std::vector<std::pair<int, int>>(simulation->num_grid_cells, std::make_pair(0, 0));

    //uniform grid
    simulation->uniform_gird_cells = std::vector<std::vector<int>> (simulation->num_grid_cells, std::vector<int>{});

    //counting sort
    simulation->counts = std::vector<int>(simulation->num_grid_cells + 1, 0);
    simulation->counting_sort_sorted_indices = std::vector<int>(simulation->num_particles, 0);

    print_resume(&simulation->bullet_physics_simulation);

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

    simulate_bullet(&simulation->bullet_physics_simulation, dt);

    dt = glm::clamp(dt, 0.001f, 0.01f);
    simulation->time_step = dt;

    int n = simulation->num_particles;
    int X = simulation->domainX;
    int Y = simulation->domainY;
    int Z = simulation->domainZ;

    //std::vector<glm::vec3>& positions = simulation->positions;
    //std::vector<glm::vec3>& positions_star = simulation->positions_star;
    std::vector<glm::vec3>& velocities = simulation->velocities;
    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;
    
    float kernelRadius = simulation->kernelRadius;

    //integration
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] += simulation->gravity * simulation->mass * dt;
        simulation->positions_star[i] = simulation->positions[i] + velocities[i] * dt; //prediction
    }

    //find_neighbors_counting_sort(simulation);
    find_neighbors_uniform_grid(simulation);
    //find_neighbors_brute_force(simulation);

    //solve pressure
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        float densitiy = 0.0;
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[neighbors[i][j]];
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
            glm::vec3 temp = simulation->positions_star[i] - simulation->positions_star[neighbors[i][j]];
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
    
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        //equation 13 (applying pressure force correction)
        //pressures_forces[i] = glm::vec3(0.0);
        glm::vec3 pressure_force = glm::vec3(0.0);
        for (int j = 0; j < neighbors[i].size(); j++) {
            glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[neighbors[i][j]];
            pressure_force += (lambdas[i] + lambdas[j] + s_coor(simulation, glm::length(ij))) * simulation->gradW(simulation, ij);
        }

        pressure_force /= simulation->rest_density;
        
        //update prediction
        simulation->positions_star[i] += pressure_force;

        //update velocity
        simulation->positions_star[i].x = resolve_collision(simulation->positions_star[i].x, simulation->particleRadius, X - simulation->particleRadius);
        simulation->positions_star[i].y = resolve_collision(simulation->positions_star[i].y, simulation->particleRadius, Y - simulation->particleRadius);
        simulation->positions_star[i].z = resolve_collision(simulation->positions_star[i].z, simulation->particleRadius, Z - simulation->particleRadius);

        velocities[i] = (simulation->positions_star[i] - simulation->positions[i]) / simulation->time_step;
        simulation->positions[i] = simulation->positions_star[i];

    }

}

};