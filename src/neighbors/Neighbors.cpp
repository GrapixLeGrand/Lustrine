
#include "Neighbors.hpp"
#include "Utils.hpp"
#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include <utility>
#include "Sorting.hpp"

namespace Lustrine {
/*
bool check_index(int i, int min, int max) {
    return (i >= min && i < max);
}*/


void clear_neighbors(Simulation* simulation) {
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        simulation->neighbors[i].clear();
    }
}

void find_neighbors_brute_force(Simulation* simulation) {

    assert(false); std::cout << "this function does not work\n" << std::endl;
    /*
    clear_neighbors(simulation);

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        glm::vec3& self = simulation->positions_star[i];
        for (int j = i + 1; j < simulation->num_particles; j++) {
            glm::vec3& other = simulation->positions_star[j];
            glm::vec3 tmp = self - other;
            if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                simulation->neighbors[i].push_back(j);
                if (j < simulation->ptr_sand_end) {
                    simulation->neighbors[j].push_back(i);
                }
            }
        } //for 2
    } // for 1
    */
}


void sort_particles_by_cell_id(Simulation* simulation, std::vector<std::pair<int, int>>& particle_cell_id) {
    std::sort(particle_cell_id.begin(), particle_cell_id.end(), 
        [](const auto& a, const auto& b)
        { 
            return a.second < b.second; 
        }
    );

    //const int lower_bound = simulation
    //for ()
}

void counting_sort_v2(Simulation* simulation) {

    memset(simulation->counting_sort_arrays->counts, 0, (simulation->num_grid_cells + 1) * sizeof(int));

    for (int i = 0; i < simulation->num_sand_particles; i++) {
        simulation->counting_sort_arrays->counts[simulation->sand_particle_cell_id[i].second]++;
    }

    for (int i = 1; i < simulation->num_grid_cells ; i++) {
        simulation->counting_sort_arrays->counts[i] += simulation->counting_sort_arrays->counts[i - 1];
    }

    for (int i = simulation->num_sand_particles - 1; i >= 0; i--) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->counting_sort_arrays->counts[cell_id]--;
        simulation->sand_particle_cell_id[simulation->counting_sort_arrays->counts[cell_id]].first = i;
        simulation->sand_particle_cell_id[simulation->counting_sort_arrays->counts[cell_id]].second = cell_id;
        //simulation->positions_star_copy[simulation->counts[cell_id]] = simulation->positions_star[i];
    }

}

void find_neighbors_uniform_grid_v3(Simulation* simulation) {

    //clear_neighbors(simulation);

    if (!simulation->computed_static_particles) {
        for (int i = simulation->ptr_solid_start; i < simulation->ptr_solid_end; i++) {
            int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
            glm::vec3& position = simulation->positions_star[i];
            simulation->uniform_grid_cells_static_saved[cell_id].push_back(i);
        }
        simulation->computed_static_particles = true;
    }

    simulation->uniform_gird_cells = simulation->uniform_grid_cells_static_saved;

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->sand_particle_cell_id[i].first = i;
        simulation->sand_particle_cell_id[i].second = cell_id;
    }

    counting_sort_v2(simulation);
    
    /*std::sort(simulation->sand_particle_cell_id.begin(), simulation->sand_particle_cell_id.end(), 
        [](const auto& a, const auto& b)
        { 
            return a.second < b.second;
        }
    );*/

    simulation->uniform_gird_cells = simulation->uniform_grid_cells_static_saved;

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        int new_index = simulation->sand_particle_cell_id[i].first;
        int cell_index = simulation->sand_particle_cell_id[i].second;
        simulation->position_neighbor_tmp[i] = simulation->positions[new_index];
        simulation->position_star_neighbor_tmp[i] = simulation->positions_star[new_index];
        simulation->velocity_tmp[i] = simulation->velocities[new_index];
        simulation->uniform_gird_cells[cell_index].push_back(i);
    }
    
    memcpy(simulation->positions, simulation->position_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->positions_star, simulation->position_star_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->velocities, simulation->velocity_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    //simulation->velocities = simulation->velocity_tmp;
    
}

void find_neighbors_uniform_grid_v2(Simulation* simulation) {

    clear_neighbors(simulation);

    if (!simulation->computed_static_particles) {
        for (int i = simulation->ptr_solid_start; i < simulation->ptr_solid_end; i++) {
            int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
            glm::vec3& position = simulation->positions_star[i];
            simulation->uniform_grid_cells_static_saved[cell_id].push_back(i);
        }
        simulation->computed_static_particles = true;
    }

    simulation->uniform_gird_cells = simulation->uniform_grid_cells_static_saved;

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        //simulation->uniform_gird_cells[cell_id].push_back(i);
        simulation->sand_particle_cell_id[i].first = i;
        simulation->sand_particle_cell_id[i].second = cell_id;
    }

    //counting_sort_v2(simulation);
    std::sort(simulation->sand_particle_cell_id.begin(), simulation->sand_particle_cell_id.end(), 
        [](const auto& a, const auto& b)
        { 
            return a.second < b.second;
        }
    );

    simulation->uniform_gird_cells = simulation->uniform_grid_cells_static_saved;

    glm::vec3 tmp {0, 0, 0};
    glm::vec3 tmp1 {0, 0, 0};
    glm::vec3 tmp2 {0, 0, 0};
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        int new_index = simulation->sand_particle_cell_id[i].first;
        int cell_index = simulation->sand_particle_cell_id[i].second;
        simulation->position_neighbor_tmp[i] = simulation->positions[new_index];
        simulation->position_star_neighbor_tmp[i] = simulation->positions_star[new_index];
        simulation->velocity_tmp[i] = simulation->velocities[new_index];
        simulation->uniform_gird_cells[cell_index].push_back(i);

    }

    memcpy(simulation->positions, simulation->position_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->positions_star, simulation->position_star_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->velocities, simulation->velocity_tmp, simulation->num_sand_particles * sizeof(glm::vec3));

    for (int yy = 0; yy < simulation->gridY; yy++) {
        for (int xx = 0; xx < simulation->gridX; xx++) {
            for (int zz = 0; zz < simulation->gridZ; zz++) {

                int cell_id = 
                    yy * simulation->gridX * simulation->gridZ +
                    xx * simulation->gridZ + 
                    zz;

                std::vector<int>& current_cell_indices = simulation->uniform_gird_cells[cell_id];

                if (current_cell_indices.empty() == true) {
                    continue;
                }

                for (int y = -1; y <= 1; y++) {
                    for (int x = -1; x <= 1; x++) {
                        for (int z = -1; z <= 1; z++) {
                            
                            if (
                                check_index(xx + x, 0, simulation->gridX) == false ||
                                check_index(yy + y, 0, simulation->gridY) == false ||
                                check_index(zz + z, 0, simulation->gridZ) == false
                            )
                            {
                                continue;
                            }

                            int neighbor_cell_id =
                                (yy + y) * simulation->gridX * simulation->gridZ +
                                (xx + x) * simulation->gridZ + 
                                (zz + z);

                            std::vector<int>& neighbor_cell_indices = simulation->uniform_gird_cells[neighbor_cell_id];

                            if (neighbor_cell_indices.empty() == true) {
                                continue;
                            }

                            for (int i = 0; i < current_cell_indices.size(); i++) {
                                const int current_index = current_cell_indices[i];
                                if (current_index >= simulation->ptr_solid_start) {
                                    continue;
                                }
                                const glm::vec3& self = simulation->positions_star[current_index];
                                for (int j = 0; j < neighbor_cell_indices.size(); j++) {
                                    const int neighbor_index = neighbor_cell_indices[j];
                                    const glm::vec3& other = simulation->positions_star[neighbor_index];
                                    const glm::vec3 tmp = self - other;
                                    if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                                        simulation->neighbors[current_index].push_back(neighbor_index);
                                    }
                                }
                            } //end distance check

                        }
                    }
                } //end neighbor cells checking 


            }
        }
    } //end grid checking
    
}


        /*tmp = simulation->positions_star[new_index];
        simulation->positions_star[new_index] = simulation->positions_star[i]; 
        simulation->positions_star[i] = tmp;

        tmp = simulation->positions[new_index];
        simulation->positions[new_index] = simulation->positions[i]; 
        simulation->positions[i] = tmp;

        tmp = simulation->velocities[new_index];
        simulation->velocities[new_index] = simulation->velocities[i]; 
        simulation->velocities[i] = tmp;*/

void find_neighbors_uniform_grid_v1(Simulation* simulation) {
    
    clear_neighbors(simulation);

    if (!simulation->computed_static_particles) {
        for (int i = simulation->ptr_solid_start; i < simulation->ptr_solid_end; i++) {
            int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
            simulation->uniform_grid_cells_static_saved[cell_id].push_back(i);
        }
        simulation->computed_static_particles = true;
    }
    
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->counting_sort_arrays->particles_unsorted_indices[i] = cell_id;
    }

    Sorting::counting_sort(
        simulation->counting_sort_arrays->counts,
        simulation->counting_sort_arrays->particles_unsorted_indices,
        simulation->num_sand_particles,
        simulation->num_grid_cells,
        simulation->counting_sort_arrays->particles_sorted_indices
    );

    memcpy(simulation->position_neighbor_tmp, simulation->positions, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->position_star_neighbor_tmp, simulation->positions_star, simulation->num_sand_particles * sizeof(glm::vec3));
    memcpy(simulation->velocity_tmp, simulation->velocities, simulation->num_sand_particles * sizeof(glm::vec3));
    simulation->uniform_gird_cells = simulation->uniform_grid_cells_static_saved;

    int* sorted_indices = simulation->counting_sort_arrays->particles_sorted_indices;
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        simulation->positions[sorted_indices[i]] = simulation->position_neighbor_tmp[i];
        simulation->positions_star[sorted_indices[i]] = simulation->position_star_neighbor_tmp[i];
        simulation->velocities[sorted_indices[i]] = simulation->velocity_tmp[i];
        int cell_id = get_cell_id(simulation, simulation->positions_star[sorted_indices[i]]);
        simulation->uniform_gird_cells[cell_id].push_back(sorted_indices[i]);

        /*
        int new_index = simulation->sand_particle_cell_id[i].first;
        int cell_index = simulation->sand_particle_cell_id[i].second;
        simulation->position_neighbor_tmp[i] = simulation->positions[new_index];
        simulation->position_star_neighbor_tmp[i] = simulation->positions_star[new_index];
        simulation->velocity_tmp[i] = simulation->velocities[new_index];
        simulation->uniform_gird_cells[cell_index].push_back(i);
        */
    }

    //memcpy(simulation->positions, simulation->position_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    //memcpy(simulation->positions_star, simulation->position_star_neighbor_tmp, simulation->num_sand_particles * sizeof(glm::vec3));
    //memcpy(simulation->velocities, simulation->velocity_tmp, simulation->num_sand_particles * sizeof(glm::vec3));


    for (int yy = 0; yy < simulation->gridY; yy++) {
        for (int xx = 0; xx < simulation->gridX; xx++) {
            for (int zz = 0; zz < simulation->gridZ; zz++) {

                int cell_id =
                    yy * simulation->gridX * simulation->gridZ +
                    xx * simulation->gridZ +
                    zz;

                std::vector<int>& current_cell_indices = simulation->uniform_gird_cells[cell_id];

                if (current_cell_indices.empty() == true) {
                    continue;
                }

                for (int y = -1; y <= 1; y++) {
                    for (int x = -1; x <= 1; x++) {
                        for (int z = -1; z <= 1; z++) {

                            if (
                                check_index(xx + x, 0, simulation->gridX) == false ||
                                check_index(yy + y, 0, simulation->gridY) == false ||
                                check_index(zz + z, 0, simulation->gridZ) == false
                                )
                            {
                                continue;
                            }

                            int neighbor_cell_id =
                                (yy + y) * simulation->gridX * simulation->gridZ +
                                (xx + x) * simulation->gridZ +
                                (zz + z);

                            std::vector<int>& neighbor_cell_indices = simulation->uniform_gird_cells[neighbor_cell_id];

                            if (neighbor_cell_indices.empty() == true) {
                                continue;
                            }

                            for (int i = 0; i < current_cell_indices.size(); i++) {
                                const int current_index = current_cell_indices[i];
                                if (current_index >= simulation->ptr_solid_start) {
                                    continue;
                                }
                                const glm::vec3& self = simulation->positions_star[current_index];
                                for (int j = 0; j < neighbor_cell_indices.size(); j++) {
                                    const int neighbor_index = neighbor_cell_indices[j];
                                    const glm::vec3& other = simulation->positions_star[neighbor_index];
                                    const glm::vec3 tmp = self - other;
                                    if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                                        simulation->neighbors[current_index].push_back(neighbor_index);
                                    }
                                }
                            } //end distance check

                        }
                    }
                } //end neighbor cells checking 


            }
        }
    } //end grid checking
}

void find_neighbors_uniform_grid(Simulation* simulation) {

    clear_neighbors(simulation);

    for (int i = 0; i < simulation->num_grid_cells; i++) {
        simulation->uniform_gird_cells[i].clear();
    }

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        glm::vec3& position = simulation->positions_star[i];
        simulation->uniform_gird_cells[cell_id].push_back(i);
    }

    for (int i = simulation->ptr_solid_start; i < simulation->ptr_solid_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        glm::vec3& position = simulation->positions_star[i];
        simulation->uniform_gird_cells[cell_id].push_back(i);
    }

    for (int yy = 0; yy < simulation->gridY; yy++) {
        for (int xx = 0; xx < simulation->gridX; xx++) {
            for (int zz = 0; zz < simulation->gridZ; zz++) {

                int cell_id = 
                    yy * simulation->gridX * simulation->gridZ +
                    xx * simulation->gridZ + 
                    zz;

                std::vector<int>& current_cell_indices = simulation->uniform_gird_cells[cell_id];

                if (current_cell_indices.empty() == true) {
                    continue;
                }

                for (int y = -1; y <= 1; y++) {
                    for (int x = -1; x <= 1; x++) {
                        for (int z = -1; z <= 1; z++) {
                            
                            if (
                                check_index(xx + x, 0, simulation->gridX) == false ||
                                check_index(yy + y, 0, simulation->gridY) == false ||
                                check_index(zz + z, 0, simulation->gridZ) == false
                            )
                            {
                                continue;
                            }

                            int neighbor_cell_id =
                                (yy + y) * simulation->gridX * simulation->gridZ +
                                (xx + x) * simulation->gridZ + 
                                (zz + z);

                            std::vector<int>& neighbor_cell_indices = simulation->uniform_gird_cells[neighbor_cell_id];

                            if (neighbor_cell_indices.empty() == true) {
                                continue;
                            }

                            for (int i = 0; i < current_cell_indices.size(); i++) {
                                const int current_index = current_cell_indices[i];
                                if (current_index >= simulation->ptr_solid_start) {
                                    continue;
                                }
                                const glm::vec3& self = simulation->positions_star[current_index];
                                for (int j = 0; j < neighbor_cell_indices.size(); j++) {
                                    const int neighbor_index = neighbor_cell_indices[j];
                                    const glm::vec3& other = simulation->positions_star[neighbor_index];
                                    const glm::vec3 tmp = self - other;
                                    if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                                        simulation->neighbors[current_index].push_back(neighbor_index);
                                    }
                                }
                            } //end distance check

                        }
                    }
                } //end neighbor cells checking 


            }
        }
    } //end grid checking

}


void find_neighbors_counting_sort(Simulation* simulation) {

    /*
    assign_particles_to_cells(simulation);
    std::sort(simulation->particle_cell_index_to_index.begin(), simulation->particle_cell_index_to_index.end(), 
        [] (const auto& p1, const auto& p2)
        {
            return p1.first < p2.first;
        }
    );*/

    /*
    auto start = std::chrono::high_resolution_clock::now();

    simulation->positions_star_copy = simulation->positions_star;
    counting_sort(simulation);
    
    //std::fill(simulation->cell_indices.begin(), simulation->cell_indices.end(), std::make_pair(0, 0));
    
    simulation->cell_indices[0].first = 0;
    simulation->cell_indices[simulation->num_particles - 1].second = simulation->num_particles;
    int current_cell_index = 0;

    for (int i = 0; i < simulation->num_particles; i++) {
        auto& sorted_particle = simulation->particle_cell_index_to_index[i];
        int sorted_cell_index = sorted_particle.first;
        //int unsorted_particle_index = sorted_particle.second;
        //simulation->positions_star_copy[i] = simulation->positions_star[unsorted_particle_index];

        if (current_cell_index != sorted_cell_index) {
            simulation->cell_indices[current_cell_index].second = i - 1;
            simulation->cell_indices[sorted_cell_index].first = i;
        }
    }

    clear_neighbors(simulation);

    for (int yy = 0; yy < simulation->gridY; yy++) {
        for (int xx = 0; xx < simulation->gridX; xx++) {
            for (int zz = 0; zz < simulation->gridZ; zz++) {
                
                int current_cell_id = 
                    yy * simulation->gridX * simulation->gridZ +
                    xx * simulation->gridZ +
                    zz;

                const auto& range = simulation->cell_indices[current_cell_id];
                int lower = range.first;
                int upper = range.second;

                if (lower == upper) {
                    continue;
                }

                for (int y = -1; y <= 1; y++) {
                    for (int x = -1; x <= 1; x++) {
                        for (int z = -1; z <= 1; z++) {
                            
                            if (
                                check_index(xx + x, 0, simulation->gridX) == false ||
                                check_index(yy + y, 0, simulation->gridY) == false ||
                                check_index(zz + z, 0, simulation->gridZ) == false
                            )
                            {
                                continue;
                            }

                            int neighbor_cell_id = 
                                (yy + y) * simulation->gridX * simulation->gridZ +
                                (xx + x) * simulation->gridZ +
                                (zz + z);

                            const auto& neighbor_range = simulation->cell_indices[neighbor_cell_id];
                            int neighbor_lower = neighbor_range.first;
                            int neighbor_upper = neighbor_range.second;

                            for (int i = lower; i < upper; i++) {
                                const glm::vec3& self = simulation->positions_star_copy[i];
                                int current_index = simulation->particle_cell_index_to_index[i].second;
                                for (int j = neighbor_lower; j < neighbor_upper; j++) {
                                    if (i == j) continue;
                                    const glm::vec3& other = simulation->positions_star_copy[j];
                                    const glm::vec3 tmp = self - other;
                                    if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                                        int neighbor_index = simulation->particle_cell_index_to_index[j].second;
                                        simulation->neighbors[current_index].push_back(neighbor_index);
                                    }
                                }
                            }
                            //end distance check
                        }
                    }
                }
                //end grid cells
            }
        }
    }

    */

    /*
    for (int i = 0; i < simulation->num_particles; i++) {
        const glm::vec3& self = simulation->positions_star_copy[i];
        //glm::vec3 indices = get_cell_id_comp(simulation, self, i);

        int sorted_cell_index = simulation->particle_cell_index_to_index[i].first;
        int unsorted_particle_index = simulation->particle_cell_index_to_index[i].second;
        
        int slice = sorted_cell_index % (simulation->gridZ * simulation->gridY);
        int xx = (sorted_cell_index - slice) / (simulation->gridZ * simulation->gridY);
        int zz = slice % simulation->gridZ;
        int yy = (slice - zz) / simulation->gridZ;

        //glm::vec3 indices = simulation->positions_star[unsorted_particle_index];
        //indices /= simulation->cell_size;
        //indices = glm::clamp(indices, glm::vec3(0.0),  glm::vec3(simulation->gridX, simulation->gridY, simulation->gridZ));
        //int xx = (int) indices.x;
        //int zz = (int) indices.y;
        //int yy = (int) indices.z;

        for (int x = -1; x <= 1; x++) {
            for (int y = -1; y <= 1; y++) {
                for (int z = -1; z <= 1; z++) {
                    
                    
                    if (
                        check_index(xx + x, 0, simulation->gridX) == false ||
                        check_index(yy + y, 0, simulation->gridY) == false ||
                        check_index(zz + z, 0, simulation->gridZ) == false
                    )
                    {
                        continue;
                    }

                    int neighor_cell_id = 
                        //sorted_cell_index + 
                        (x + xx) * simulation->gridY * simulation->gridZ +
                        (y + yy) * simulation->gridZ +
                        (z + zz);

                    if (neighor_cell_id >= simulation->num_grid_cells || neighor_cell_id < 0) {
                        continue;
                    }

                    neighor_cell_id = glm::clamp(neighor_cell_id, 0, simulation->num_grid_cells - 1);

                    const auto& range = simulation->cell_indices[neighor_cell_id];
                    int lower = range.first;
                    int upper = range.second;
                    for (int k = lower; k < upper; k++) {
                        if (k == i) continue;
                        int unsorted_other_particle_index = simulation->particle_cell_index_to_index[k].second;
                        const glm::vec3& other = simulation->positions_star_copy[k];
                        const glm::vec3 tmp = self - other;
                        if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                            simulation->neighbors[unsorted_particle_index].push_back(unsorted_other_particle_index);
                        }
                    }

                }
            }
        }

    }*/
}

}; //namespace Lustrine