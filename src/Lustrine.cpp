#include "Lustrine.hpp"

#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include "profiling/Profiling.hpp"
#include "Simulate.hpp"

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
    const std::vector<Grid>& grids_solid_arg
) {
    init_simulation(
        parameters, 
        simulation, 
        grids_sand_arg, 
        grids_solid_arg,
        1
    );
}

void init_simulation(
    const SimulationParameters* parameters, 
    Simulation* simulation, 
    const std::vector<Grid>& grids_sand_arg, 
    const std::vector<Grid>& grids_solid_arg,
    int subdivision
) {

    Profiling::init_profiling();
    Bullet::init_bullet(&simulation->bullet_physics_simulation);
    
    simulation->simulate_fun = simulate_sand_v1;

    std::vector<glm::vec3> grids_sand_positions_arg;
    std::vector<glm::vec3> grids_solid_positions_arg;

    grids_sand_positions_arg.reserve(grids_sand_arg.size());
    for (int i = 0; i < grids_sand_arg.size(); i++) {
        grids_sand_positions_arg.push_back(grids_sand_arg[i].position);
    }

    grids_solid_positions_arg.reserve(grids_solid_arg.size());
    for (int i = 0; i < grids_solid_arg.size(); i++) {
        grids_solid_positions_arg.push_back(grids_solid_arg[i].position);
    }

    simulation->domainX = parameters->X;
    simulation->domainY = parameters->Y;
    simulation->domainZ = parameters->Z;

    int first_guess_allocation = simulation->domainX * simulation->domainY * simulation->domainZ * subdivision * subdivision * subdivision;
    assert(first_guess_allocation > 0);
    int initial_allocated_particles_num = 0;

    for (int i = 0; i < grids_sand_arg.size(); i++) {
        initial_allocated_particles_num += grids_sand_arg[i].num_occupied_grid_cells;
    }

    simulation->num_sand_particles = initial_allocated_particles_num;

    for (int i = 0; i < grids_solid_arg.size(); i++) {
        initial_allocated_particles_num += grids_solid_arg[i].num_occupied_grid_cells * subdivision * subdivision * subdivision;
    }

    simulation->num_solid_particles = initial_allocated_particles_num - simulation->num_sand_particles;

    simulation->total_allocated = (size_t) std::max(initial_allocated_particles_num, first_guess_allocation);
    simulation->leftover_allocated = simulation->total_allocated - initial_allocated_particles_num;

    simulation->positions = new glm::vec3[simulation->total_allocated];
    simulation->positions_star = new glm::vec3[simulation->total_allocated];
    simulation->colors = new glm::vec4[simulation->total_allocated];
    simulation->positions_tmp = new glm::vec3[simulation->total_allocated];

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
            init_chunk_from_grid(parameters, &chunk, &simulation->grids_sand[i], SAND);
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
    }// end sand

    {//solids
        
        simulation->ptr_solid_ordered_start = simulation->total_allocated - 1; //simulation->total_allocated - simulation->num_solid_particles; //simulation->total_allocated - 1; //head is going down so list is in order
        simulation->ptr_solid_ordered_end = simulation->total_allocated - 1;

        simulation->grids_solid = grids_solid_arg;
        simulation->grids_initial_positions_solid = grids_solid_positions_arg;
        simulation->chunks_solid.reserve(simulation->grids_solid.size());
        simulation->solid_grid_to_body = std::vector<std::pair<int, int>>(simulation->grids_solid.size(), std::make_pair(-1, -1));
        simulation->grids_solid_chunk_ptrs = std::vector<std::pair<int, int>>(simulation->grids_solid.size(), std::make_pair(-1, -1));

        for (int i = 0; i < simulation->grids_solid.size(); i++) {
            assert(simulation->grids_solid[i].type == SOLID);
            Chunk chunk;
            init_chunk_from_grid_subdivision(parameters, &chunk, &simulation->grids_solid[i], SOLID, subdivision);
            //init_chunk_from_grid(parameters, &chunk, &simulation->grids_solid[i], SOLID);
            simulation->chunks_solid.push_back(chunk);
            simulation->grids_solid_chunk_ptrs[i].second = simulation->ptr_solid_ordered_end + 1;

            // add boxes in bullet
            Chunk chunk_boxes;
            init_chunk_from_grid_unit_length(parameters, &chunk_boxes, &simulation->grids_solid[i], SOLID);
            for (int j = 0; j < chunk_boxes.num_particles; j++) {
                assert(simulation->grids_solid[i].dynamic_solid == false);
                int bullet_body_index = Bullet::add_box(&simulation->bullet_physics_simulation, chunk_boxes.positions[j], simulation->grids_solid[i].dynamic_solid);
            }
            for (int i = 0; i < 10; ++i) {
                std::cout << chunk.positions[i].x << " " << chunk.positions[i].y << " " << chunk.positions[i].z << std::endl;
            }

            for (int j = 0; j < chunk.num_particles; j++) {
                
                assert(simulation->grids_solid[i].dynamic_solid == false);

                simulation->positions[simulation->ptr_solid_ordered_end] = chunk.positions[j];
                simulation->positions_star[simulation->ptr_solid_ordered_end] = chunk.positions[j];
                if (chunk.has_one_color_per_particles == true) {
                    simulation->colors[simulation->ptr_solid_ordered_end] = chunk.colors[j];
                } else {
                    simulation->colors[simulation->ptr_solid_ordered_end] = chunk.color;
                }
                simulation->ptr_solid_ordered_end--;
            }
            simulation->grids_solid_chunk_ptrs[i].first = simulation->ptr_solid_ordered_end;
        }

        simulation->ptr_solid_start = simulation->ptr_solid_ordered_end + 1;
        simulation->ptr_solid_end = simulation->ptr_solid_ordered_start + 1; //start was inclusive

        if (simulation->ptr_solid_ordered_start == simulation->ptr_solid_ordered_end) { //not sure
            simulation->ptr_solid_start = simulation->ptr_solid_ordered_start;
            simulation->ptr_solid_end = simulation->ptr_solid_ordered_end;
        }

        //TODO rethink these asserts...
        //assert(simulation->ptr_solid_start <= simulation->ptr_solid_ordered_start);
        //assert(simulation->ptr_solid_end >= simulation->ptr_solid_ordered_start + 1);
        //assert(simulation->ptr_solid_start >= simulation->ptr_solid_ordered_end + 1);

        assert(simulation->ptr_solid_ordered_start == simulation->total_allocated - 1);
        assert(simulation->ptr_solid_ordered_end == simulation->ptr_solid_ordered_start - simulation->num_solid_particles);

        simulation->positions_solid = simulation->positions + simulation->total_allocated - simulation->num_solid_particles;
        simulation->colors_solid = simulation->colors + simulation->total_allocated - simulation->num_solid_particles;
    
    }// end solids

    simulation->velocities = std::vector<glm::vec3>(simulation->num_sand_particles, {0.0f, 0.0f, 0.0f});
    simulation->lambdas = std::vector<float>(simulation->num_sand_particles, 0.0f);
    simulation->neighbors = std::vector<std::vector<int>>(simulation->num_sand_particles, std::vector<int>{});

    simulation->W = cubic_kernel;
    simulation->gradW = cubic_kernel_grad;

    //sorting neighbor strategy with grid
    simulation->particleRadius = parameters->particleRadius;
    simulation->particleDiameter = parameters->particleDiameter;
    simulation->kernelRadius = 3.1f * parameters->particleRadius;
    simulation->cell_size = 1.0f * simulation->kernelRadius;

    //for the kernel
    float h3 = std::pow(simulation->kernelRadius, 3);
    simulation->cubic_kernel_k = 8.0f / (pi * h3);
    simulation->cubic_kernel_l = 48.0f / (pi * h3);

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

    Bullet::allocate_particles_colliders(&simulation->bullet_physics_simulation, simulation->bullet_physics_simulation.num_particles_allocated, simulation->particleRadius);
    Bullet::bind_foreign_sand_positions(&simulation->bullet_physics_simulation, simulation->positions);
    Bullet::print_resume(&simulation->bullet_physics_simulation);

}

void clean_simulation(Simulation* simulation) {
    clean_bullet(&simulation->bullet_physics_simulation);
    delete[] simulation->positions;
    delete[] simulation->positions_star;
    delete[] simulation->colors;
    delete[] simulation->positions_tmp;
}

void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, glm::vec3 position, glm::vec4 color, MaterialType type) {

    grid->type = type;
    grid->num_grid_cells = X * Y * Z;
    grid->num_occupied_grid_cells = grid->num_grid_cells;
    grid->X = X;
    grid->Y = Y;
    grid->Z = Z;
    grid->has_one_color_per_cell = false;

    grid->cells = std::vector<bool>(grid->num_grid_cells, true);
    grid->color = color;
    
    grid->position = position;

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;

    grid->sparse_solid = true;
    grid->dynamic_solid = false;

}

void init_chunk_from_grid(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, MaterialType type) {
    
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

                    particle_position += grid->position;

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


void init_chunk_from_grid_subdivision(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, MaterialType type, int subdivision) {

    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells * subdivision * subdivision * subdivision; // 225
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, { 0, 0, 0 });
    if (chunk->has_one_color_per_particles) {
        chunk->colors = std::vector<glm::vec4>(chunk->num_particles, { 0, 0, 0, 1 });
    }
    else {
        chunk->color = grid->color;
    }

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;
    glm::vec3 offset = glm::vec3(radius, radius, radius);
    int X = grid->X * subdivision; int Y = grid->Y * subdivision; int Z = grid->Z * subdivision;
    int counter = 0; //because chunk is disorganized we cannot use grid indices

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                int cell_index = (x / subdivision) * grid->Y * grid->Z + (y / subdivision) * grid->Z + (z / subdivision);
                if (grid->cells[cell_index] == true) { //a particle present at the current grid cell
                    glm::vec3& particle_position = chunk->positions[counter];

                    particle_position.x = x * diameter;
                    particle_position.y = y * diameter;
                    particle_position.z = z * diameter;

                    particle_position += grid->position;

                    if (chunk->has_one_color_per_particles == true) {
                        const glm::vec4& grid_color = grid->colors[cell_index];
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

void init_chunk_from_grid_unit_length(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, MaterialType type) {

    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells; // 225
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, { 0, 0, 0 });
    if (chunk->has_one_color_per_particles) {
        chunk->colors = std::vector<glm::vec4>(chunk->num_particles, { 0, 0, 0, 1 });
    }
    else {
        chunk->color = grid->color;
    }

    const float diameter = 1.0f;
    const float radius = 0.5f;
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

                    particle_position += grid->position;

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

//void init_chunk_from_grid_subdivision(const SimulationParameters* parameters, Chunk* chunk, const Grid* grid, MaterialType type, int subdivision) {


void simulate(Simulation* simulation, float dt) {
    Profiling::start_counter(0);
    simulation->simulate_fun(simulation, dt);
    Profiling::stop_counter(0);
}




};