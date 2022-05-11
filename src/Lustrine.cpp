#include "Lustrine.hpp"

#include <iostream>
#include <chrono>
#include <set>
#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include "profiling/Profiling.hpp"
#include "Simulate.hpp"
#include "neighbors/Utils.hpp"

namespace Lustrine {

constexpr double pi = 3.14159265358979323846;


void add_physics_static_grid(Simulation* simulation, const Grid* grid) {

    for (int x = 0; x < grid->X; x++) {
        for (int y = 0; y < grid->Y; y++) {
            for (int z = 0; z < grid->Z; z++) {
                int index = grid->Y * grid->Z * x + grid->Z * y + z;
                if (grid->cells[index]) {
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

    simulation->simulate_fun = simulate_sand;

    simulation->subdivision = subdivision;
    simulation->parameters_copy = *parameters;

    simulation->source = new ParticleSource();
    simulation->source->num_sources = 0;
    simulation->sink = new ParticleSink();
    simulation->sink->num_sinks = 0;
    simulation->sink->temp_removal.reserve(100000);//safe amount to avoid resize
    simulation->wind_system = new WindSystem();
    simulation->wind_system->direction = glm::vec3(0, 0, -1);
    simulation->wind_system->magnitude = 5.0f;

    simulation->bullet_physics_simulation.particle_radius = parameters->particleRadius;

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

    std::align_val_t simd_vector_align{64};
    simulation->positions = new (simd_vector_align)glm::vec3[simulation->total_allocated];
    simulation->positions_star = new (simd_vector_align) glm::vec3[simulation->total_allocated];
    simulation->colors = new (simd_vector_align) glm::vec4[simulation->total_allocated];
    simulation->positions_tmp = new (simd_vector_align) glm::vec3[simulation->total_allocated];

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
            init_chunk_from_grid(&chunk, &simulation->grids_sand[i], SAND, parameters->particleDiameter, 1, false);
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
            init_chunk_from_grid(&chunk, &simulation->grids_solid[i], SOLID, parameters->particleDiameter, subdivision, true);
            simulation->chunks_solid.push_back(chunk);
            simulation->grids_solid_chunk_ptrs[i].second = simulation->ptr_solid_ordered_end + 1;

            // add boxes in bullet
            Chunk chunk_boxes;
            init_chunk_from_grid(&chunk_boxes, &simulation->grids_solid[i], SOLID, 1.0f, 1, false);
            for (int j = 0; j < chunk_boxes.num_particles; j++) {
                assert(simulation->grids_solid[i].dynamic_solid == false);
                int bullet_body_index = Bullet::add_box(&simulation->bullet_physics_simulation, chunk_boxes.positions[j], simulation->grids_solid[i].dynamic_solid);
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

    simulation->num_remaining_sand_particles = simulation->ptr_solid_ordered_end - simulation->ptr_sand_end;//TODO check this

    int total_possible_num_sand_particles = simulation->ptr_solid_ordered_start;

    simulation->velocities = new (simd_vector_align) glm::vec3[total_possible_num_sand_particles];
    memset(simulation->velocities, 0, total_possible_num_sand_particles);
    //simulation->velocities = std::vector<glm::vec3>(simulation->num_sand_particles, {0.0f, 0.0f, 0.0f});
    simulation->lambdas = std::vector<float>(total_possible_num_sand_particles, 0.0f);
    simulation->neighbors = std::vector<std::vector<int>>(total_possible_num_sand_particles, std::vector<int>{});

    simulation->W = cubic_kernel;
    simulation->gradW = cubic_kernel_grad;

    //sorting neighbor strategy with grid
    simulation->particleRadius = parameters->particleRadius;
    simulation->particleDiameter = parameters->particleDiameter;
    simulation->kernelRadius = 3.1f * parameters->particleRadius;//be careful with these parameters they affect both stabiliy and performance! Do not change without checking performance and stability
    simulation->cell_size = 1.0f * simulation->kernelRadius;

    //for the kernel
    float h3 = std::pow(simulation->kernelRadius, 3);
    simulation->cubic_kernel_k = 8.0f / (pi * h3);
    simulation->cubic_kernel_l = 48.0f / (pi * h3);

    simulation->gridX = (int) (simulation->domainX / simulation->cell_size) + 1; //plus 1 because its a full range
    simulation->gridY = (int) (simulation->domainY / simulation->cell_size) + 1;
    simulation->gridZ = (int) (simulation->domainZ / simulation->cell_size) + 1;

    //uniform grid
    simulation->num_grid_cells = simulation->gridX * simulation->gridY * simulation->gridZ;
    simulation->uniform_gird_cells = std::vector<std::vector<int>> (simulation->num_grid_cells, std::vector<int>{});

    //perf test
    simulation->uniform_grid_cells_static_saved = std::vector<std::vector<int>> (simulation->num_grid_cells, std::vector<int>{});
    simulation->sand_particle_cell_id = std::vector<std::pair<int, int>>(total_possible_num_sand_particles, std::make_pair(0, 0));
    simulation->position_neighbor_tmp = new (simd_vector_align) glm::vec3[total_possible_num_sand_particles];
    simulation->position_star_neighbor_tmp = new (simd_vector_align) glm::vec3[total_possible_num_sand_particles];
    simulation->velocity_tmp = new (simd_vector_align) glm::vec3[total_possible_num_sand_particles]; //std::vector<glm::vec3>(simulation->num_sand_particles, {0, 0, 0});
    //simulation->colors_neighbor_tmp = new glm::vec4[simulation->num_sand_particles];
    //perf test

    {
        //used in the counting sort
        simulation->counting_sort_arrays = new CountingSortArrays();
        simulation->counting_sort_arrays->counts = new (simd_vector_align) int[simulation->num_grid_cells + 1];
        simulation->counting_sort_arrays->particles_unsorted_indices = new (simd_vector_align) int[total_possible_num_sand_particles];
        simulation->counting_sort_arrays->particles_sorted_indices = new (simd_vector_align) int[total_possible_num_sand_particles];
    }

    Bullet::allocate_particles_colliders(&simulation->bullet_physics_simulation, simulation->bullet_physics_simulation.num_particles_allocated, simulation->particleRadius);
    Bullet::bind_foreign_sand_positions(&simulation->bullet_physics_simulation, simulation->positions);
    Bullet::print_resume(&simulation->bullet_physics_simulation);

    std::cout << "registered " << simulation->num_sand_particles << " sand particles and " << simulation->num_solid_particles << " solid particles\n";
    std::cout << "registered (at start at lease)" << simulation->bullet_physics_simulation.ptr_bounding_box_end - simulation->bullet_physics_simulation.ptr_bounding_box_start << " boxes for colliding with the player\n";
    std::cout << "sand start " << simulation->ptr_sand_start << "\n";
    std::cout << "sand end " << simulation->ptr_sand_end << "\n";
    std::cout << "solid start " << simulation->ptr_solid_start << "\n";
    std::cout << "solid end " << simulation->ptr_solid_end << "\n";

}

void clean_simulation(Simulation* simulation) {
    clean_bullet(&simulation->bullet_physics_simulation);

    std::align_val_t simd_vector_align{64};
    ::operator delete[] (simulation->positions, simd_vector_align);
    ::operator delete[] (simulation->positions_star, simd_vector_align);
    ::operator delete[] (simulation->colors, simd_vector_align);
    ::operator delete[] (simulation->positions_tmp, simd_vector_align);
    ::operator delete[] (simulation->velocities, simd_vector_align);

    ::operator delete[] (simulation->position_neighbor_tmp, simd_vector_align);
    ::operator delete[] (simulation->position_star_neighbor_tmp, simd_vector_align);
    ::operator delete[] (simulation->velocity_tmp, simd_vector_align);

    delete simulation->source;
    delete simulation->sink;
    delete simulation->wind_system;

    ::operator delete[] (simulation->counting_sort_arrays->counts, simd_vector_align);
    ::operator delete[] (simulation->counting_sort_arrays->particles_sorted_indices, simd_vector_align);
    ::operator delete[] (simulation->counting_sort_arrays->particles_unsorted_indices, simd_vector_align);
    delete simulation->counting_sort_arrays;
}

void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, glm::vec3 position, glm::vec4 color, MaterialType type) {

    grid->type = type;
    grid->num_grid_cells = X * Y * Z;
    grid->num_occupied_grid_cells = grid->num_grid_cells;
    grid->X = X;
    grid->Y = Y;
    grid->Z = Z;
    grid->has_one_color_per_cell = false;

    grid->cells = std::vector<int>(grid->num_grid_cells, true);
    grid->color = color;
    
    grid->position = position;

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;

    grid->sparse_solid = true;
    grid->dynamic_solid = false;

}


void init_grid_box_random(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, glm::vec3 position, glm::vec4 color, MaterialType type, float probability) {

    grid->type = type;
    grid->num_grid_cells = X * Y * Z;
    grid->num_occupied_grid_cells = 0;
    grid->X = X;
    grid->Y = Y;
    grid->Z = Z;
    grid->has_one_color_per_cell = false;

    grid->cells = std::vector<int>(grid->num_grid_cells, true);

    for (int y = 0; y < Y; y++) {
        for (int x = 0; x < X; x++) {
            for (int z = 0; z < Z; z++) {
                float p = ((float) rand()) / RAND_MAX;
                grid->cells[y * X * Z + x * Z + z] = (p <= probability);
                grid->num_occupied_grid_cells++;
            }
        }
    }

    grid->color = color;
    
    grid->position = position;

    const float diameter = parameters->particleDiameter;
    const float radius = parameters->particleRadius;

    grid->sparse_solid = true;
    grid->dynamic_solid = false;

}

void init_chunk_from_grid(Chunk* chunk, const Grid* grid, MaterialType type, float cell_size, int subdivision, bool mask_out) {

    chunk->type = type;
    chunk->num_particles = grid->num_occupied_grid_cells * subdivision * subdivision * subdivision;
    chunk->has_one_color_per_particles = grid->has_one_color_per_cell;

    chunk->positions = std::vector<glm::vec3>(chunk->num_particles, { 0, 0, 0 });
    if (chunk->has_one_color_per_particles) {
        chunk->colors = std::vector<glm::vec4>(chunk->num_particles, { 0, 0, 0, 1 });
    }
    else {
        chunk->color = grid->color;
    }

    glm::vec3 offset = glm::vec3(cell_size * 0.5f, cell_size * 0.5f, cell_size * 0.5f);
    int X = grid->X * subdivision; int Y = grid->Y * subdivision; int Z = grid->Z * subdivision;
    int counter = 0; //because chunk is disorganized we cannot use grid indices

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                int cell_index = (x / subdivision) * grid->Y * grid->Z + (y / subdivision) * grid->Z + (z / subdivision);
                if (grid->cells[cell_index] && ((!mask_out) || (mask_out && grid->cells[cell_index]!=1))) { //a particle present at the current grid cell
                    glm::vec3& particle_position = chunk->positions[counter];

                    particle_position.x = x * cell_size;
                    particle_position.y = y * cell_size;
                    particle_position.z = z * cell_size;

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


float total_time = 0.0f;

int ll = 0;
float kk = 0.0f;

void simulate(Simulation* simulation, float dt) {

/*
    if (ll % 2 == 0) {
        for (int i = 0; i < 100000; i++) {
            kk += std::sqrt(kk);
        }
        std::cout << "doo the work " << dt << "\n";
    } else {
        std::cout << "nope " << dt << "\n"; 
    }

    ll++;
*/
float dt_clamped = 0.016f; //glm::clamp(dt, 0.001f, 0.016f);
    for (int i = 0; i < simulation->source->num_sources; i++) {
        simulation->source->timers[i] += dt_clamped;
    }

    for (int i = 0; i < simulation->source->num_sources; i++) {
        if (simulation->source->spawned[i] < simulation->source->capacities[i] && simulation->source->timers[i] >= simulation->source->frequencies[i]) {
            Chunk& pattern = simulation->source->patterns[i];
            if (simulation->num_remaining_sand_particles < pattern.num_particles && simulation->source->source_state[i]) {
                continue;
            }
            
            float freq = simulation->source->frequencies[i];
            float t_last = total_time;- simulation->source->timers[i];
            float offset = freq - fmodf(t_last, freq);  
            int num_to_spawn = std::floor((simulation->source->timers[i] + offset) / freq);
            
            float particle_diameter = 2.0f * simulation->particleRadius;
            float speed = freq > 0.0f ? 1.0f * particle_diameter / freq : 1.0f;
            
            std::cout << "total time: " << total_time << ", t_last:" << t_last << ", spawed " << num_to_spawn << " offset:" << offset << "\n";
            glm::vec3& direction = simulation->source->directions[i];
            glm::vec3 offset_freq = direction * offset; //fractional offset
            
            for (int j = 0; j < 1; j++) {
                
                glm::vec3 current_offset = j * particle_diameter * direction + offset_freq;
                for (int k = 0; k < pattern.num_particles; k++) {
                    simulation->positions[k + simulation->ptr_sand_end] = pattern.positions[k] + current_offset;
                    simulation->velocities[k + simulation->ptr_sand_end] = direction * speed + simulation->gravity * (j * freq);
                    simulation->colors[k + simulation->ptr_sand_end] = pattern.color;
                }
                
            }

            //memcpy(simulation->positions + simulation->ptr_sand_end, pattern.positions.data(), pattern.num_particles * sizeof(glm::vec3));
            //memcpy(simulation->positions_star + simulation->ptr_sand_end, pattern.positions.data(), pattern.num_particles * sizeof(glm::vec3));
            //memset(simulation->colors + simulation->ptr_sand_end, simulation->spawner.patterns[i].color);

            int end_ptr_tmp = simulation->ptr_sand_end;
            simulation->ptr_sand_end += pattern.num_particles;
            simulation->num_remaining_sand_particles -= pattern.num_particles;
            simulation->num_sand_particles += pattern.num_particles;

            /*
            glm::vec3& direction = simulation->source->directions[i];
            float factor = simulation->source->frequencies[i] > 0.0f ? 1.5f * (simulation->particleRadius * 2.0f) / simulation->source->frequencies[i] : 1.0f;
            for (int j = end_ptr_tmp; j < simulation->ptr_sand_end; j++) {
                simulation->velocities[j] = direction * factor;
                simulation->colors[j] = simulation->source->patterns[i].color;
            }

            simulation->velocities[end_ptr_tmp + pattern.num_particles / 2] += 0.001f * glm::vec3(1.0); //small pertuabation to avoid particles stacking
            simulation->source->timers[i] = 0.0f;
            simulation->source->spawned[i] += pattern.num_particles;
            */
           simulation->velocities[end_ptr_tmp + pattern.num_particles / 2] += 0.001f * glm::vec3(1.0); //small pertuabation to avoid particles stacking
            simulation->source->timers[i] = 0.0f;
            simulation->source->spawned[i] += pattern.num_particles;
        }
    }

    Profiling::start_counter(0);
    simulation->simulate_fun(simulation, dt);
    Profiling::stop_counter(0);

    for (int i = 0; i < simulation->sink->num_sinks; i++) {

        if (simulation->sink->state[i] && simulation->sink->timers[i] >= simulation->sink->frequencies[i]) {
            int current_despawned = 0;
            const int sink_cells_size = simulation->sink->sink_cells[i].size();
            for (int j = 0; j < sink_cells_size; j++) {
                int cell_id = simulation->sink->sink_cells[i][j];
                int cell_size = simulation->uniform_gird_cells[cell_id].size();
                for (int k = 0; k < cell_size; k++) {
                    int evicted_id = simulation->uniform_gird_cells[cell_id][k];
                    if (evicted_id < simulation->ptr_sand_end)
                        simulation->sink->temp_removal.push_back(evicted_id);
                }
                current_despawned += cell_size;
            }
        }

    }

    for (int i = 0; i < simulation->sink->temp_removal.size(); i++) {
        int evicted = simulation->sink->temp_removal[i];
        simulation->positions[evicted] = simulation->positions[simulation->ptr_sand_end - 1];
        simulation->positions_star[evicted] = simulation->positions_star[simulation->ptr_sand_end - 1];
        simulation->velocities[evicted] = simulation->velocities[simulation->ptr_sand_end - 1];
        simulation->ptr_sand_end--;
        simulation->num_sand_particles--;
        simulation->num_remaining_sand_particles++;
    }

    simulation->sink->temp_removal.clear();

    /*
            std::vector<int>& cell = simulation->uniform_gird_cells[i];
        for (int j = 0; j < cell.size(); j++) {
            int evicted = cell[j];
            if (cell[j] < simulation->ptr_sand_end) {

                simulation->positions[evicted] = simulation->positions[simulation->ptr_sand_end - 1];
                simulation->positions_star[evicted] = simulation->positions_star[simulation->ptr_sand_end - 1];
                simulation->velocities[evicted] = simulation->velocities[simulation->ptr_sand_end - 1];
                simulation->ptr_sand_end--;
                simulation->num_sand_particles--;
                simulation->num_remaining_sand_particles++;

            }
        }*/

    

    for (int i = 0; i < simulation->sink->num_sinks; i++) {
        simulation->sink->timers[i] += dt;
    }

    total_time += dt_clamped;
}


int query_cell_num_particles(Simulation* simulation, glm::vec3 min_pos, glm::vec3 max_pos, bool include_solid) {

    int counter = 0;
    int min_x = 0, min_y = 0, min_z = 0;
    int max_x = 0, max_y = 0, max_z = 0;

    CellIndices min_indices = get_cell_indices(simulation, min_pos);
    min_x = (int) min_indices.x;
    min_y = (int) min_indices.y;
    min_z = (int) min_indices.z;

    CellIndices max_indices = get_cell_indices(simulation, max_pos);

    max_x = (int) max_indices.x;
    max_y = (int) max_indices.y;
    max_z = (int) max_indices.z;

    if (include_solid) {
        for (int y = min_y; y <= max_y; y++) {
            for (int x = min_x; x <= max_x; x++) {
                for (int z = min_z; z <= max_z; z++) {
                    int cell_id =
                        y * simulation->gridX * simulation->gridZ +
                        x * simulation->gridZ +
                        z;
                    counter += simulation->uniform_gird_cells[cell_id].size();
                }
            }
        }
    } else {
        for (int y = min_y; y <= max_y; y++) {
            for (int x = min_x; x <= max_x; x++) {
                for (int z = min_z; z <= max_z; z++) {
                    int cell_id =
                        y * simulation->gridX * simulation->gridZ +
                        x * simulation->gridZ +
                        z;
                    
                    const int size = simulation->uniform_gird_cells[cell_id].size();
                    for (int i = 0; i < size; i++) {
                        if (simulation->uniform_gird_cells[cell_id][i] < simulation->ptr_sand_end) {
                            counter++;
                        }
                    }
                }
            }
        }
    }

    return counter;
}


int add_particle_source(Simulation* simulation, const Grid* pattern, glm::vec3 direction, float freq, int capacity) {

    Chunk chunk;
    init_chunk_from_grid(&chunk, pattern, Lustrine::MaterialType::SAND, simulation->particleDiameter, simulation->subdivision, false);
    int index = simulation->source->num_sources;
    simulation->source->num_sources++;
    simulation->source->patterns.push_back(chunk);
    direction /= glm::length(direction);
    simulation->source->directions.push_back(direction);
    simulation->source->timers.push_back(0.0);
    simulation->source->frequencies.push_back(freq);
    simulation->source->source_state.push_back(true);
    if (capacity < 0) {
        capacity = std::numeric_limits<int>::max();
    }
    simulation->source->capacities.push_back(capacity);
    simulation->source->spawned.push_back(0);
    return index;
}

int add_particle_sink(Simulation* simulation, glm::vec3 min_pos, glm::vec3 max_pos, float frequency) {
    
    int min_x = 0, min_y = 0, min_z = 0;
    int max_x = 0, max_y = 0, max_z = 0;

    CellIndices min_indices = get_cell_indices(simulation, min_pos);
    min_x = (int) min_indices.x;
    min_y = (int) min_indices.y;
    min_z = (int) min_indices.z;

    CellIndices max_indices = get_cell_indices(simulation, max_pos);
    max_x = (int) max_indices.x;
    max_y = (int) max_indices.y;
    max_z = (int) max_indices.z;

    int index = simulation->sink->num_sinks;
    simulation->sink->sink_cells.push_back(std::vector<int>{});
    simulation->sink->despawned.push_back(0);
    simulation->sink->frequencies.push_back(frequency);
    simulation->sink->state.push_back(true);
    simulation->sink->timers.push_back(0.0f);

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            for (int z = min_z; z <= max_z; z++) {
                int cell_id =
                    y * simulation->gridX * simulation->gridZ +
                    x * simulation->gridZ +
                    z;
                simulation->sink->sink_cells[index].push_back(cell_id);
            }
        }
    }

    simulation->sink->num_sinks++;
    std::sort(simulation->sink->sink_cells[index].begin(), simulation->sink->sink_cells[index].end());
    return index;
}



int add_particle_sink(Simulation* simulation, const Grid* pattern, float frequency) {

    int index = simulation->sink->num_sinks;
    simulation->sink->sink_cells.push_back(std::vector<int>{});
    simulation->sink->despawned.push_back(0);
    simulation->sink->frequencies.push_back(frequency);
    simulation->sink->state.push_back(true);
    simulation->sink->timers.push_back(0.0f);

    std::set<int> indices;

    std::cout << "sink add with grid not implemented yet!" << std::endl;

    /*
    for (int x = 0; x < pattern->X; x++) {
        for (int y = 0; y < pattern->Y; y++) {
            for (int z = 0; z < pattern->Z; z++) {
                int index = 
                pattern->
            }
        }
    }*/
    return 0;
}

void set_source_state(Simulation* simulation, int index, bool state) {
    simulation->source->source_state[index] = state;
}

void set_sink_state(Simulation* simulation, int index, bool state) {
    simulation->sink->state[index] = state;
}

int get_source_spawned(Simulation* simulation, int index) {
    return simulation->source->spawned[index];
}

int get_sink_despawned(Simulation* simulation, int index) {
    return simulation->sink->despawned[index];
}

/*
void remove_particle_sink(Simulation* simulation, int index, glm::vec3 min_pos, glm::vec3 max_pos) {
    int min_x = 0, min_y = 0, min_z = 0;
    int max_x = 0, max_y = 0, max_z = 0;

    CellIndices min_indices = get_cell_indices(simulation, min_pos);
    min_x = (int) min_indices.x;
    min_y = (int) min_indices.y;
    min_z = (int) min_indices.z;

    CellIndices max_indices = get_cell_indices(simulation, max_pos);

    max_x = (int) max_indices.x;
    max_y = (int) max_indices.y;
    max_z = (int) max_indices.z;

    std::vector<int> temp;
    std::vector<int> sink_cells_diff;
    sink_cells_diff.reserve(simulation->sink->sink_cells.size());
    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            for (int z = min_z; z <= max_z; z++) {
                int cell_id =
                    y * simulation->gridX * simulation->gridZ +
                    x * simulation->gridZ +
                    z;
                temp.push_back(cell_id);
            }
        }
    }

    std::sort(temp.begin(), temp.end());
    std::set_difference(
        simulation->sink->sink_cells.begin(),
        simulation->sink->sink_cells.end(), 
        temp.begin(),
        temp.end(),
        std::back_inserter(sink_cells_diff)
    );

    simulation->sink->sink_cells = sink_cells_diff;

}*/

void update_wind_system(WindSystem* wind_system, float dt) {

    

}

};