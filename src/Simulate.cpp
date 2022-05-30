#include "Simulate.hpp"
#include "profiling/Profiling.hpp"
#include "neighbors/Neighbors.hpp"

namespace Lustrine {

    float s_coor(const Simulation* simulation, float rl) {
        return -simulation->s_corr_k * std::pow(simulation->W(simulation, rl) / simulation->W(simulation, simulation->s_corr_dq), simulation->s_corr_n);
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


    void simulate_fluid(Simulation* simulation, float dt) {

        simulate_bullet(&simulation->bullet_physics_simulation, dt, simulation->ptr_sand_start, simulation->ptr_sand_end);

        dt = glm::clamp(dt, 0.001f, 0.01f);
        simulation->time_step = dt;

        int n = simulation->num_particles;
        int X = simulation->domainX;
        int Y = simulation->domainY;
        int Z = simulation->domainZ;

        //std::vector<glm::vec3>& positions = simulation->positions;
        //std::vector<glm::vec3>& positions_star = simulation->positions_star;
        glm::vec3* velocities = simulation->velocities;
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
                glm::vec3 neighbor_grad = -(simulation->mass / simulation->rest_density) * simulation->gradW(simulation, temp);
                constraint_gradient_sum += glm::dot(neighbor_grad, neighbor_grad);
                grad_current_p -= neighbor_grad;
            }

            constraint_gradient_sum += glm::dot(grad_current_p, grad_current_p);

            lambdas[i] = 0.0;
            if (constraint_gradient_sum > 0.0) {
                lambdas[i] = -constraint_i / (constraint_gradient_sum + simulation->relaxation_epsilon);
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



    inline glm::vec3 solve_boundary_collision_constraint(glm::vec3 n, glm::vec3 p0, glm::vec3 p, float d) {
        float C = glm::dot(n, (p - p0)) - d;
        if (C >= 0) {
            return glm::vec3(0.0);
        }

        // https://matthias-research.github.io/pages/publications/posBasedDyn.pdf Eq(9)
        glm::vec3 dC = n;
        float s = C / glm::dot(dC, dC);
        glm::vec3 dp = -s * dC;

        return dp;
    }

    // require x >= 0
    inline float avoid0(float x) {
        return x + 1e-9f;
    }

    float attract_kernel(float x) {
        return 1.0f;
    }


    float blow_kernel(float x, float kernel_radius) {
        float alpha = 0.75f;
        if (x > kernel_radius) {
            return 0.0f;
        }
        else {
            float value = 1.0f - x / kernel_radius;

            return value;
        }
        //return x;
    }

void simulate_sand(Simulation* simulation, float dt) {
    Bullet::simulate_bullet(&simulation->bullet_physics_simulation, dt, simulation->ptr_sand_start, simulation->ptr_sand_end);

    float collision_coeff = 0.8f;
    float boundary_collision_coeff = 0.9f;
    float friction_coeff = 0.7f;
    float mu_s = 0.95f;
    float mu_k = 0.8f;
    
    //dt = 0.016;
    //dt = glm::clamp(dt, 0.001f, 0.016f); //TEMPORARY
    //dt = (1.0f / 30.0f);
    simulation->time_step = dt;

    int n = simulation->num_particles;
    float X = simulation->domainX;
    float Y = simulation->domainY;
    float Z = simulation->domainZ;

    glm::vec3* positions = simulation->positions;
    glm::vec3* positions_star = simulation->positions_star;
    glm::vec3* velocities = simulation->velocities;


    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;

    float kernelRadius = simulation->kernelRadius;

    static bool prev_attract_flag = false;

    //integration
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        float w = 1.0f / simulation->mass;
        velocities[i] += simulation->gravity * dt;
        if (prev_attract_flag && !simulation->attract_flag) {
            simulation->attracted[i] = false;
        }
        if (simulation->attract_flag) {
            if (glm::length(simulation->bullet_physics_simulation.player_position - positions[i]) < simulation->attract_radius) {
                simulation->attracted[i] = true;
            }
            if (simulation->attracted[i]) {
                glm::vec3 particle_to_attraction_origin = (simulation->bullet_physics_simulation.player_position + glm::vec3(0.0f, 1.5f, 0.0f)) - positions[i]; // above player's head
                velocities[i] += glm::normalize(particle_to_attraction_origin)* attract_kernel(glm::length(particle_to_attraction_origin)) * simulation->attract_coeff * simulation->particleRadius * dt * w;
            }
        }
        if (simulation->blow_flag) {
            glm::vec3 particle_to_blow_origin = simulation->bullet_physics_simulation.player_position - positions[i];
            if (glm::length(particle_to_blow_origin) < simulation->blow_radius) {
                velocities[i] += -glm::normalize(particle_to_blow_origin) * blow_kernel(glm::length(particle_to_blow_origin), simulation->blow_radius) * simulation->blow_coeff * simulation->particleRadius * w;
            }
        }

        positions_star[i] = positions[i] + velocities[i] * dt; 
        // clamp particle positions to be inside boundaries
        positions_star[i] = glm::clamp(positions_star[i], glm::vec3(simulation->particleRadius), glm::vec3(X, Y, Z) - glm::vec3(simulation->particleRadius));
    }
    //glm::vec3 *positions_tmp = new glm::vec3[simulation->total_allocated];
    find_neighbors_uniform_grid_v1(simulation);
    //memcpy(simulation->positions_tmp, positions_star, sizeof(glm::vec3) * simulation->total_allocated);
    if (simulation->first_iteration) {
        memcpy(simulation->positions_tmp, positions_star, sizeof(glm::vec3) * simulation->total_allocated);
        simulation->first_iteration = false;
    } else {
        memcpy(simulation->positions_tmp, positions_star, sizeof(glm::vec3) * simulation->num_sand_particles);
    }

        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    for (int substep = 0; substep < 4; ++substep) {

        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
            glm::vec3 deltap = glm::vec3(0, 0, 0);
            for (int j : neighbors[i]) {
                if (i == j) continue;
                glm::vec3 pi = simulation->positions_tmp[i];
                glm::vec3 pj = simulation->positions_tmp[j];

                glm::vec3 ij = pi - pj;
                if (glm::length(ij) == 0.0f) {
                    ij = glm::vec3(0.0f, 0.00001f, 0.0f);
                }
                float len = glm::length(ij);
                if (len > simulation->particleDiameter) continue;

                if (j >= simulation->ptr_sand_start && j < simulation->ptr_sand_end)
                {
                    float mass = simulation->mass;
                    float neighbor_mass = simulation->mass;
                    glm::vec3 tmp = collision_coeff * (neighbor_mass) / (mass + neighbor_mass) * (len - simulation->particleDiameter) * ij / len;

                    deltap -= tmp;

                    float d = glm::length(tmp);
                    glm::vec3 xidelta = (pi - tmp) - positions[i];
                    glm::vec3 xjdelta = (pj + tmp) - positions[j];
                    glm::vec3 normalX = (pi - tmp) - (pj + tmp);
                    normalX = glm::normalize(normalX);
                    glm::vec3 relativeX = xidelta - xjdelta;

                    glm::vec3 xtan = (relativeX - dot(relativeX, normalX) * normalX);
                    if ((d * mu_s) > avoid0(glm::length(xtan))) {
                        deltap -= friction_coeff * xtan;
                    }
                    else {
                        float ratio = (mu_k * d / avoid0(glm::length(xtan))) < 1 ? (mu_k * d / avoid0(glm::length(xtan))) : 1;
                        deltap -= friction_coeff * xtan * ratio;
                    }
                }
                else
                {
                    glm::vec3 tmp = collision_coeff * (len - simulation->particleDiameter) * ij / len;
                    deltap -= tmp;

                    float d = length(tmp);
                    glm::vec3 xidelta = (pi - tmp) - positions[i];
                    glm::vec3 xjdelta = glm::vec3(0, 0, 0);
                    glm::vec3 normalX = (pi - tmp) - pj;
                    normalX = glm::normalize(normalX);
                    glm::vec3 relativeX = xidelta - xjdelta;

                    glm::vec3 xtan = (relativeX - dot(relativeX, normalX) * normalX);
                    if ((d * mu_s) > avoid0(glm::length(xtan))) {
                        deltap -= friction_coeff * xtan;
                    }
                    else {
                        float ratio = (mu_k * d / avoid0(glm::length(xtan))) < 1 ? (mu_k * d / avoid0(glm::length(xtan))) : 1;
                        deltap -= friction_coeff * xtan * ratio;
                    }
                }
            }
            positions_star[i] = simulation->positions_tmp[i] + deltap;


            //here
        }


        // particle-bounadry collision
        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
            glm::vec3 p = simulation->positions_tmp[i];

            //glm::vec3 dp = glm::vec3(0.0);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(1, 0, 0), glm::vec3(simulation->particleRadius, 0, 0), p, simulation->particleRadius);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(-1, 0, 0), glm::vec3(X- simulation->particleRadius, 0, 0), p, simulation->particleRadius);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 1, 0), glm::vec3(0, simulation->particleRadius, 0), p, simulation->particleRadius);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, -1, 0), glm::vec3(0, Y- simulation->particleRadius, 0), p, simulation->particleRadius);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 0, 1), glm::vec3(0, 0, simulation->particleRadius), p, simulation->particleRadius);
            //dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 0, -1), glm::vec3(0, 0, Z- simulation->particleRadius), p, simulation->particleRadius);
            //positions_star[i] += dp;
            positions_star[i] = glm::clamp(positions_star[i], glm::vec3(simulation->particleRadius), glm::vec3(X, Y, Z) - glm::vec3(simulation->particleRadius));
        }

        memcpy(simulation->positions_tmp + simulation->ptr_sand_start, positions_star + simulation->ptr_sand_start, sizeof(glm::vec3) * simulation->num_sand_particles);
    }



    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] = (positions_star[i] - positions[i]) / simulation->time_step;
        positions[i] = positions_star[i];
    }
    Profiling::stop_counter(6);


    prev_attract_flag = simulation->attract_flag;

}


void simulate_sand_v3(Simulation* simulation, float dt) {
    Bullet::simulate_bullet(&simulation->bullet_physics_simulation, dt, simulation->ptr_sand_start, simulation->ptr_sand_end);

    float collision_coeff = 0.8f;
    float boundary_collision_coeff = 0.9f;
    float friction_coeff = 0.5f;
    float mu_s = 0.8f;
    float mu_k = 0.8f;

    dt = glm::clamp(dt, 0.001f, 0.016f);
    //dt = (1.0f / 30.0f);
    simulation->time_step = dt;

    int n = simulation->num_particles;
    float X = simulation->domainX;
    float Y = simulation->domainY;
    float Z = simulation->domainZ;

    glm::vec3* positions = simulation->positions;
    glm::vec3* positions_star = simulation->positions_star;
    glm::vec3* velocities = simulation->velocities;


    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;

    float kernelRadius = simulation->kernelRadius;

    //integration
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        float w = 1.0f / simulation->mass;
        velocities[i] += simulation->gravity * dt;
        if (simulation->attract_flag) {
            glm::vec3 p_to_cp = simulation->bullet_physics_simulation.player_position - positions[i];
            //        float force_magnitude = 10 * std::min(1.0, 1.0/ glm::dot(p_to_gs, p_to_gs));
            velocities[i] += glm::normalize(p_to_cp) * simulation->particleRadius * 1000.0f * dt * w;
        }
        if (simulation->blow_flag) {
            glm::vec3 p_to_cp = simulation->bullet_physics_simulation.player_position - positions[i];
            //        float force_magnitude = 10 * std::min(1.0, 1.0/ glm::dot(p_to_gs, p_to_gs));
            velocities[i] += -glm::normalize(p_to_cp) * simulation->particleRadius * 200.0f * w;
        }

        positions_star[i] = positions[i] + velocities[i] * dt; // update both
    }
    //glm::vec3 *positions_tmp = new glm::vec3[simulation->total_allocated];
    
    
    find_neighbors_uniform_grid_v1(simulation);
    memcpy(simulation->positions_tmp, positions_star, sizeof(glm::vec3) * simulation->total_allocated);

        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {

        for (int substep = 0; substep < 4; ++substep) {
            glm::vec3 deltap = glm::vec3(0, 0, 0);
            for (int j : neighbors[i]) {
                if (i == j) continue;
                glm::vec3 pi = simulation->positions_tmp[i];
                glm::vec3 pj = simulation->positions_tmp[j];

                glm::vec3 ij = pi - pj;
                float len = glm::length(ij);
                if (len > simulation->particleDiameter) continue;

                if (j >= simulation->ptr_sand_start && j < simulation->ptr_sand_end)
                {
                    float mass = simulation->mass;
                    float neighbor_mass = simulation->mass;
                    glm::vec3 tmp = collision_coeff * (neighbor_mass) / (mass + neighbor_mass) * (len - simulation->particleDiameter) * ij / len;

                    deltap -= tmp;

                    float d = glm::length(tmp);
                    glm::vec3 xidelta = (pi - tmp) - positions[i];
                    glm::vec3 xjdelta = (pj + tmp) - positions[j];
                    glm::vec3 normalX = (pi - tmp) - (pj + tmp);
                    normalX = glm::normalize(normalX);
                    glm::vec3 relativeX = xidelta - xjdelta;

                    glm::vec3 xtan = (relativeX - dot(relativeX, normalX) * normalX);
                    if ((d * mu_s) > avoid0(glm::length(xtan))) {
                        deltap -= friction_coeff * xtan;
                    }
                    else {
                        float ratio = (mu_k * d / avoid0(glm::length(xtan))) < 1 ? (mu_k * d / avoid0(glm::length(xtan))) : 1;
                        deltap -= friction_coeff * xtan * ratio;
                    }
                }
                else
                {
                    glm::vec3 tmp = collision_coeff * (len - simulation->particleDiameter) * ij / len;
                    deltap -= tmp;

                    float d = length(tmp);
                    glm::vec3 xidelta = (pi - tmp) - positions[i];
                    glm::vec3 xjdelta = glm::vec3(0, 0, 0);
                    glm::vec3 normalX = (pi - tmp) - pj;
                    normalX = glm::normalize(normalX);
                    glm::vec3 relativeX = xidelta - xjdelta;

                    glm::vec3 xtan = (relativeX - dot(relativeX, normalX) * normalX);
                    if ((d * mu_s) > avoid0(glm::length(xtan))) {
                        deltap -= friction_coeff * xtan;
                    }
                    else {
                        float ratio = (mu_k * d / avoid0(glm::length(xtan))) < 1 ? (mu_k * d / avoid0(glm::length(xtan))) : 1;
                        deltap -= friction_coeff * xtan * ratio;
                    }
                }
            }
            positions_star[i] = simulation->positions_tmp[i] + deltap;

        }
    }

    
        // particle-bounadry collision
        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
            glm::vec3 p = simulation->positions_tmp[i];
            float r = simulation->particleRadius;

            glm::vec3 dp = glm::vec3(0.0);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(1, 0, 0), glm::vec3(0, 0, 0), p, r);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(-1, 0, 0), glm::vec3(X, 0, 0), p, r);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 1, 0), glm::vec3(0, 0, 0), p, r);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, -1, 0), glm::vec3(0, Y, 0), p, r);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), p, r);
            dp += boundary_collision_coeff * solve_boundary_collision_constraint(glm::vec3(0, 0, -1), glm::vec3(0, 0, Z), p, r);
            positions_star[i] += dp;
        }
    

        memcpy(simulation->positions_tmp + simulation->ptr_sand_start, positions_star + simulation->ptr_sand_start, sizeof(glm::vec3) * simulation->num_sand_particles);
    



    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] = (positions_star[i] - positions[i]) / simulation->time_step;
        positions[i] = positions_star[i];
    }
    Profiling::stop_counter(6);

    //TODO: this perturbation is for test, remove this in the future
    /*static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.01, 0.01, 0.01);
        perturbation = false;
    }*/

}

void simulate_sand_v1(Simulation* simulation, float dt) {

    Profiling::start_counter(1);
    Bullet::simulate_bullet(&simulation->bullet_physics_simulation, dt, simulation->ptr_sand_start, simulation->ptr_sand_end);
    Profiling::stop_counter(1);

    Profiling::start_counter(2);

    float collision_coeff = 0.9f;
    float boundary_collision_coeff = 0.9f;
    float friction_coeff = 0.5f;
    float mu_s = 0.5f;
    float mu_k = 0.8f;

    dt = (1.0f / 30.0f);
    simulation->time_step = dt;

    int X = simulation->domainX;
    int Y = simulation->domainY;
    int Z = simulation->domainZ;

    glm::vec3* positions = simulation->positions;
    glm::vec3* positions_star = simulation->positions_star;
    glm::vec3* velocities = simulation->velocities;

    std::vector<std::vector<int>>& neighbors = simulation->neighbors;

    //integration

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] += simulation->gravity * dt;
        positions_star[i] = positions[i] + velocities[i] * dt; // update both
    }

    //find_neighbors_counting_sort(simulation);
    Profiling::start_counter(3);
    find_neighbors_uniform_grid_v2(simulation);
    Profiling::stop_counter(3);

    float particleDiameter2 = simulation->particleDiameter * simulation->particleDiameter;
        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    float mass_inv = 1.0f / simulation->mass;
    Profiling::start_counter(4);
    const int lower_bound = simulation->ptr_sand_start;
    const int upper_bound = simulation->ptr_sand_end;

    for (int substep = 0; substep < 4; ++substep) {

        //std::vector<std::pair<int, int>> contacts;
        for (int i = lower_bound; i < upper_bound; i++) {
            for (int j : neighbors[i]) {
                if (i >= j) continue;
                glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[j];
                float len2 = glm::dot(ij, ij); //glm::length(ij);
                if (len2 < particleDiameter2) {

                    float len = glm::length(ij);

                    // TODO: can be optimized by using a single mass;
                    // i is always sand
                    float wi = mass_inv;
                    // compute j's mass
                    float wj = (j >= simulation->ptr_sand_start && j < simulation->ptr_sand_end) ? mass_inv : 0.0f;
                    float wi_wj_inv = 1.0f / (wi + wj);
                    float wi_wiPwj = wi * wi_wj_inv;
                    float wj_wiPwj = wj * wi_wj_inv;

                    glm::vec3 delta_ij = (len - simulation->particleDiameter) * ij / avoid0(len);
                    glm::vec3 delta_pi = -wi_wiPwj * delta_ij;
                    glm::vec3 delta_pj = +wj_wiPwj * delta_ij;
                    glm::vec3 positions_star_i = simulation->positions_star[i] + collision_coeff * delta_pi;
                    glm::vec3 positions_star_j = simulation->positions_star[j] + collision_coeff * delta_pj;

                    float d = simulation->particleDiameter - len;
                    glm::vec3 n = glm::normalize(positions_star_i - positions_star_j);
                    // FIXME: Should use initial positions but currently very unstable
                    // glm::vec3 delta_x_star_ij = (positions_star_i - positions[i]) - (positions_star_j - positions[j]);
                    glm::vec3 delta_x_star_ij = (positions_star_i - simulation->positions_star[i]) - (positions_star_j - simulation->positions_star[j]);
                    glm::vec3 delta_x_tangent = delta_x_star_ij - glm::dot(delta_x_star_ij, n) * n;
                    // https://mmacklin.com/uppfrta_preprint.pdf, eq(24)
                    glm::vec3 delta_pij_friction;
                    if (glm::length(delta_x_tangent) < mu_s * d) {
                        delta_pij_friction = delta_x_tangent;
                    }
                    else {
                        delta_pij_friction = delta_x_tangent * std::min(mu_k * d / avoid0(glm::length(delta_x_tangent)), 1.0f);
                    }
                    glm::vec3 delta_pi_friction = wi_wiPwj * delta_pij_friction;
                    glm::vec3 delta_pj_friction = -wj_wiPwj * delta_pij_friction;
                    positions_star_i += friction_coeff * delta_pi_friction;
                    positions_star_j += friction_coeff * delta_pi_friction;

                    positions_star[i] += positions_star_i - simulation->positions_star[i];
                    positions_star[j] += positions_star_j - simulation->positions_star[j];

                }
            }
        }


        Profiling::start_counter(5);
        // particle-bounadry collision
        for (int i = lower_bound; i < upper_bound; i++) {
            glm::vec3 p = simulation->positions_star[i];
            float r = simulation->particleRadius;

            glm::vec3 dp = glm::vec3(0.0);
            dp += solve_boundary_collision_constraint(glm::vec3(1, 0, 0), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(-1, 0, 0), glm::vec3(X, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 1, 0), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, -1, 0), glm::vec3(0, Y, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 0, -1), glm::vec3(0, 0, Z), p, r);
            positions_star[i] += boundary_collision_coeff * dp;
        }
        Profiling::stop_counter(5);
        //memcpy(simulation->positions_tmp + simulation->ptr_sand_start, positions_star + simulation->ptr_sand_start, sizeof(glm::vec3) * simulation->num_sand_particles);


    }

    Profiling::stop_counter(4);


    Profiling::start_counter(6);
    float dt_inv = 1.0f / simulation->time_step;
    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] = (positions_star[i] - positions[i]) * dt_inv;
        positions[i] = positions_star[i];
    }
    Profiling::stop_counter(6);

    Profiling::stop_counter(2);
    /*
    std::cout << "hello" << std::endl;
    static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.1, 0.1, 0.1);
        perturbation = false;
    }*/

    static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.01, 0.01, 0.01);
        perturbation = false;
    }
}


void simulate_sand_v2(Simulation* simulation, float dt) {

    Profiling::start_counter(1);
    Bullet::simulate_bullet(&simulation->bullet_physics_simulation, dt, simulation->ptr_sand_start, simulation->ptr_sand_end);
    Profiling::stop_counter(1);

    Profiling::start_counter(2);

    float collision_coeff = 0.9f;
    float boundary_collision_coeff = 0.9f;
    float friction_coeff = 0.5f;
    float mu_s = 0.5f;
    float mu_k = 0.8f;

    dt = (1.0f / 30.0f);
    simulation->time_step = dt;

    int X = simulation->domainX;
    int Y = simulation->domainY;
    int Z = simulation->domainZ;

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        simulation->velocities[i] += simulation->gravity * dt;
        simulation->positions_star[i] = simulation->positions[i] + simulation->velocities[i] * dt; // update both
    }

    //find_neighbors_counting_sort(simulation);
    Profiling::start_counter(3);
    find_neighbors_uniform_grid_v3(simulation);
    Profiling::stop_counter(3);

    float particleDiameter2 = simulation->particleDiameter * simulation->particleDiameter;
        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    float mass_inv = 1.0f / simulation->mass;
    Profiling::start_counter(4);
    const int lower_bound = simulation->ptr_sand_start;
    const int upper_bound = simulation->ptr_sand_end;

    for (int substep = 0; substep < 4; ++substep) {

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

                    int y_lower = -1; int y_upper = 1;
                    int x_lower = -1; int x_upper = 1;
                    int z_lower = -1; int z_upper = 1;
                    if (yy + y_lower < 0) y_lower = 0; if (yy + y_upper >= simulation->gridY) y_upper = 0;
                    if (xx + x_lower < 0) x_lower = 0; if (xx + x_upper >= simulation->gridX) x_upper = 0;
                    if (zz + z_lower < 0) z_lower = 0; if (zz + z_upper >= simulation->gridZ) z_upper = 0;

                    for (int y = y_lower; y <= y_upper; y++) {
                        for (int x = x_lower; x <= x_upper; x++) {
                            for (int z = z_lower; z <= z_upper; z++) {

                                int neighbor_cell_id =
                                    (yy + y) * simulation->gridX * simulation->gridZ +
                                    (xx + x) * simulation->gridZ + 
                                    (zz + z);

                                std::vector<int>& neighbor_cell_indices = simulation->uniform_gird_cells[neighbor_cell_id];

                                if (neighbor_cell_indices.empty() == true) {
                                    continue;
                                }

                                for (int ii = 0; ii < current_cell_indices.size(); ii++) {
                                    const int i = current_cell_indices[ii];
                                    if (i >= simulation->ptr_solid_start) {
                                        continue;
                                    }
                                    //const glm::vec3& self = simulation->positions_star[current_index];
                                    for (int jj = 0; jj < neighbor_cell_indices.size(); jj++) {                                        
                                        const int j = neighbor_cell_indices[jj];
                                        if (i >= j) continue;
                                        glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[j];
                                        float len2 = glm::dot(ij, ij); //glm::length(ij);
                                        if (len2 < particleDiameter2) {

                                            float len = glm::length(ij);

                                            // TODO: can be optimized by using a single mass;
                                            // i is always sand
                                            float wi = mass_inv;
                                            // compute j's mass
                                            float wj = (j >= simulation->ptr_sand_start && j < simulation->ptr_sand_end) ? mass_inv : 0.0f;
                                            float wi_wj_inv = 1.0f / (wi + wj);
                                            float wi_wiPwj = wi * wi_wj_inv;
                                            float wj_wiPwj = wj * wi_wj_inv;

                                            glm::vec3 delta_ij = (len - simulation->particleDiameter) * ij / avoid0(len);
                                            glm::vec3 delta_pi = -wi_wiPwj * delta_ij;
                                            glm::vec3 delta_pj = +wj_wiPwj * delta_ij;
                                            glm::vec3 positions_star_i = simulation->positions_star[i] + collision_coeff * delta_pi;
                                            glm::vec3 positions_star_j = simulation->positions_star[j] + collision_coeff * delta_pj;

                                            float d = simulation->particleDiameter - len;
                                            glm::vec3 n = glm::normalize(positions_star_i - positions_star_j);
                                            // FIXME: Should use initial positions but currently very unstable
                                            // glm::vec3 delta_x_star_ij = (positions_star_i - positions[i]) - (positions_star_j - positions[j]);
                                            glm::vec3 delta_x_star_ij = (positions_star_i - simulation->positions_star[i]) - (positions_star_j - simulation->positions_star[j]);
                                            glm::vec3 delta_x_tangent = delta_x_star_ij - glm::dot(delta_x_star_ij, n) * n;
                                            // https://mmacklin.com/uppfrta_preprint.pdf, eq(24)
                                            glm::vec3 delta_pij_friction;
                                            if (glm::length(delta_x_tangent) < mu_s * d) {
                                                delta_pij_friction = delta_x_tangent;
                                            }
                                            else {
                                                delta_pij_friction = delta_x_tangent * std::min(mu_k * d / avoid0(glm::length(delta_x_tangent)), 1.0f);
                                            }
                                            glm::vec3 delta_pi_friction = wi_wiPwj * delta_pij_friction;
                                            glm::vec3 delta_pj_friction = -wj_wiPwj * delta_pij_friction;
                                            positions_star_i += friction_coeff * delta_pi_friction;
                                            positions_star_j += friction_coeff * delta_pi_friction;

                                            simulation->positions_star[i] += positions_star_i - simulation->positions_star[i];
                                            simulation->positions_star[j] += positions_star_j - simulation->positions_star[j];

                                        }
                                    } // end in neighbor cell 
                                } // end in current cell 

                            }
                        }
                    } //end neighbor cells checking 


                }
            }
        } //end grid checking

        Profiling::start_counter(5);
        // particle-bounadry collision
        for (int i = lower_bound; i < upper_bound; i++) {
            glm::vec3 p = simulation->positions_star[i];
            float r = simulation->particleRadius;

            glm::vec3 dp = glm::vec3(0.0);
            dp += solve_boundary_collision_constraint(glm::vec3(1, 0, 0), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(-1, 0, 0), glm::vec3(X, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 1, 0), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, -1, 0), glm::vec3(0, Y, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 0, 1), glm::vec3(0, 0, 0), p, r);
            dp += solve_boundary_collision_constraint(glm::vec3(0, 0, -1), glm::vec3(0, 0, Z), p, r);
            simulation->positions_star[i] += boundary_collision_coeff * dp;
        }
        Profiling::stop_counter(5);
        //memcpy(simulation->positions_tmp + simulation->ptr_sand_start, positions_star + simulation->ptr_sand_start, sizeof(glm::vec3) * simulation->num_sand_particles);


    }
    Profiling::stop_counter(4);


    Profiling::start_counter(6);
    float dt_inv = 1.0f / simulation->time_step;
    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        simulation->velocities[i] = (simulation->positions_star[i] - simulation->positions[i]) * dt_inv;
        simulation->positions[i] = simulation->positions_star[i];
    }
    Profiling::stop_counter(6);

    Profiling::stop_counter(2);
    /*
    std::cout << "hello" << std::endl;
    static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.1, 0.1, 0.1);
        perturbation = false;
    }*/
    
    static bool perturbation = true;
    if (perturbation) {
        simulation->velocities[simulation->ptr_sand_start] += glm::vec3(0.01, 0.01, 0.01);
        perturbation = false;
    }
}

}