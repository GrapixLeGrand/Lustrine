#include "Lustrine.hpp"

#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>

namespace Lustrine {

constexpr double pi = 3.14159265358979323846;

void init_sim(Simulation* simulation, int particlesX, int particlesY, int particlesZ) {

    simulation->particlesX = particlesX;
    simulation->particlesY = particlesY;
    simulation->particlesZ = particlesZ;
    simulation->num_particles = particlesX * particlesY * particlesZ;

    simulation->max_neighbors = simulation->num_particles - 1;

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
}

void fill_grid(Simulation* simulation) {

    int X = simulation->particlesX;
    int Y = simulation->particlesY;
    int Z = simulation->particlesZ;

    float diameter = simulation->particleDiameter;
    float radius = simulation->particleRadius;
    glm::vec3 offset = glm::vec3(radius, radius, radius);

    for (int x = 0; x < X; x++) {
        for (int y = 0; y < Y; y++) {
            for (int z = 0; z < Z; z++) {
                
                glm::vec3& position = simulation->positions[x * Y * Z + y * Z + z];
                glm::vec4& color = simulation->colors[x * Y * Z + y * Z + z];

                position.x = x * diameter;
                position.y = y * diameter; 
                position.z = z * diameter; 

                color.r = position.x / X;
                color.g = position.y / Y;
                color.b = position.z / Z;
                color.a = 1.0f;

                position += offset;

            }
        }
    }

}

void clear_neighbors(Simulation* simulation) {
    for (int i = 0; i < simulation->num_particles; i++) {
        simulation->neighbors[i].clear();
    }
}

void find_neighbors(Simulation* simulation) {

    clear_neighbors(simulation);

    for (int i = 0; i < simulation->num_particles; i++) {
        glm::vec3& self = simulation->positions_star[i];
        for (int j = i + 1; j < simulation->num_particles; j++) {
            glm::vec3& other = simulation->positions_star[j];
            glm::vec3 tmp = self - other;
            if (glm::dot(tmp, tmp) <= simulation->kernelRadius * simulation->kernelRadius) {
                simulation->neighbors[i].push_back(j);
                simulation->neighbors[j].push_back(i);
            }
        }
    }

}


float cubic_kernel(const Simulation* simulation, float r) {
    float q = (r * simulation->kernelFactor) / simulation->kernelRadius;
    float result = 0.0;
    if (q <= 1.0) {
        if (q <= 0.5) {
            float q2 = q * q;
            float q3 = q2 * q;
            result = simulation->cubic_kernel_k * (6.0f * q3 - 6.0f * q2 + 1.0f);
        } else {
            result = simulation->cubic_kernel_k * (2.0f * std::pow(1.0f - q, 3.0f));
        }
    }
    return result;
}

float cubic_kernel(const Simulation* simulation, glm::vec3& r) {
    float tmp = glm::length(r);
    return cubic_kernel(simulation, tmp);
}

glm::vec3 cubic_kernel_grad(const Simulation* simulation, const glm::vec3& r) {
    glm::vec3 result = glm::vec3(0.0);
    float rl = glm::length(r) * simulation->kernelFactor;
    float q = rl / simulation->kernelRadius;

    if (rl > 1.0e-5 && q <= 1.0) {
        const glm::vec3 grad_q = (1.0f / (rl * simulation->kernelRadius)) * r;
        if (q <= 0.5) {
            result = simulation->cubic_kernel_l * q * (3.0f * q - 2.0f) * grad_q;
        } else {
            const float f = 1.0f - q;
            result = simulation->cubic_kernel_l * (- f * f) * grad_q;
        }
    }
    return result;
}

float poly6_kernel(const Simulation* simulation, float r) {
    float result = 0.0;
    float hf = simulation->kernelRadius * simulation->kernelFactor;
    if (r <= simulation->kernelRadius) {
        result = (315.0f / (64.0f * 3.14f * std::pow(hf, 9))) * 
            std::pow(std::pow(hf, 2) - std::pow(simulation->kernelFactor * r, 2), 3);
    }
    return result;
}

float poly6_kernel(const Simulation* simulation, glm::vec3& r) {
    return poly6_kernel(simulation, r.length());
}

glm::vec3 spiky_kernel(const Simulation* simulation, glm::vec3& r) {
    glm::vec3 result = glm::vec3(0.0);
    float rl = glm::length(r);
    if (rl > 0.0 && rl <= simulation->kernelRadius) {
        float hf = simulation->kernelRadius * simulation->kernelFactor;
        float temp = ((15.0f / (3.14f * std::pow(hf, 6))) * 
            std::pow(hf - (rl * simulation->kernelFactor), 2));
        result = (r / (rl * simulation->kernelFactor)) * temp;
    }
    return result;
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

    //glm::vec3 position = simulation->positions_star[i];
    position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
    position /= simulation->cell_size;
    //position += 1;
    /*
    int cell_id =
            simulation->gridY * simulation->gridZ + 
            ((int) position.x) * simulation->gridY * simulation->gridZ +
            simulation->gridZ + 
            ((int) position.y) * simulation->gridZ +
            1 +
            ((int) position.z);*/
    int cell_id =
            (std::floor(position.y)) * simulation->gridX * simulation->gridZ +
            (std::floor(position.x)) * simulation->gridZ +
            (std::floor(position.z));

    //cell_id = glm::clamp(cell_id, 0, simulation->num_grid_cells - 1);
    //if (cell_id < 0 || cell_id >= simulation->num_grid_cells) {
    //    std::cout << "no" << std::endl;
    //}

    return cell_id;
}

void assign_particles_to_cells(Simulation* simulation) {
    simulation->particle_cell_index_to_index.clear();
    for (int i = 0; i < simulation->num_particles; i++) {
        int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
        simulation->particle_cell_index_to_index.push_back(std::make_pair(cell_id, i));
    }
}


bool check_index(int i, int min, int max) {
    return (i >= min && i < max);
}


void find_neighbors_uniform_grid(Simulation* simulation) {

    clear_neighbors(simulation);

    for (int i = 0; i < simulation->num_grid_cells; i++) {
        simulation->uniform_gird_cells[i].clear();
    }

    for (int i = 0; i < simulation->num_particles; i++) {
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

    assign_particles_to_cells(simulation);
    std::sort(simulation->particle_cell_index_to_index.begin(), simulation->particle_cell_index_to_index.end(), 
        [] (const auto& p1, const auto& p2)
        {
            return p1.first < p2.first;
        }
    );

    simulation->positions_star_copy = simulation->positions_star;
    std::fill(simulation->cell_indices.begin(), simulation->cell_indices.end(), std::make_pair(0, 0));
    
    simulation->cell_indices[0].first = 0;
    simulation->cell_indices[simulation->num_particles - 1].second = simulation->num_particles;
    int current_cell_index = 0;

    for (int i = 0; i < simulation->num_particles; i++) {
        auto& sorted_particle = simulation->particle_cell_index_to_index[i];
        int sorted_cell_index = sorted_particle.first;
        int unsorted_particle_index = sorted_particle.second;
        simulation->positions_star_copy[i] = simulation->positions_star[unsorted_particle_index];

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

        glm::vec3 indices = simulation->positions_star[unsorted_particle_index];
        indices /= simulation->cell_size;
        indices = glm::clamp(indices, glm::vec3(0.0),  glm::vec3(simulation->gridX, simulation->gridY, simulation->gridZ));
        int xx = (int) indices.x;
        int zz = (int) indices.y;
        int yy = (int) indices.z;

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

void simulate(Simulation* simulation) {

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
    float dt = simulation->time_step;

    //integration
    for (int i = 0; i < n; i++) {
        velocities[i] += simulation->gravity * simulation->mass * dt;
        positions_star[i] = positions[i] + velocities[i] * dt; //prediction
    }

    find_neighbors_counting_sort(simulation);
    //find_neighbors_uniform_grid(simulation);
    //find_neighbors(simulation);

    //solve pressure
    for (int i = 0; i < n; i++) {

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
    
    for (int i = 0; i < n; i++) {

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

};