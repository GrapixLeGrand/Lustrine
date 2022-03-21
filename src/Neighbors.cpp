
#include "Neighbors.hpp"

#include <iostream>
#include <chrono>

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace Lustrine {

bool check_index(int i, int min, int max) {
    return (i >= min && i < max);
}

glm::vec3 get_cell_id_comp(const Simulation* simulation, glm::vec3 position, int i) {
    //glm::vec3 position = simulation->positions_star[i];
    position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
    position /= simulation->cell_size;
    position.x = (int) position.x;
    position.y = (int) position.y;
    position.z = (int) position.z;
    return position;
}

int get_cell_id(const Simulation* simulation, glm::vec3 position) {

    position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
    position /= simulation->cell_size;
    int cell_id =
            (std::floor(position.y)) * simulation->gridX * simulation->gridZ +
            (std::floor(position.x)) * simulation->gridZ +
            (std::floor(position.z));
            
    return cell_id;
}


void assign_particles_to_cells(Simulation* simulation) {
    simulation->particle_cell_index_to_index.clear();
    //std::fill(particle_cell_index_to_index.begin(), particle_cell_index_to_index.end(), std::make_pair(0, 0));
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->particle_cell_index_to_index.push_back(std::make_pair(cell_id, i));
    }

    for (int i = simulation->ptr_solid_start; i < simulation->ptr_solid_end; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->particle_cell_index_to_index.push_back(std::make_pair(cell_id, i));
    }
}


//fills the index array position star copy with position star
void counting_sort(Simulation* simulation) {

    std::fill(simulation->counts.begin(), simulation->counts.end(), 0);
    std::fill(simulation->particle_cell_index_to_index.begin(), simulation->particle_cell_index_to_index.end(), std::make_pair(0, 0));

    for (int i = 0; i < simulation->num_particles; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->counts[cell_id]++;
    }

    for (int i = 1; i < simulation->counts.size(); i++) {
        simulation->counts[i] += simulation->counts[i - 1];
    }

    for (int i = simulation->num_particles - 1; i >= 0; i--) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->counts[cell_id]--;
        simulation->particle_cell_index_to_index[simulation->counts[cell_id]].first = cell_id;
        simulation->particle_cell_index_to_index[simulation->counts[cell_id]].second = i;
        simulation->positions_star_copy[simulation->counts[cell_id]] = simulation->positions_star[i];
    }

}

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