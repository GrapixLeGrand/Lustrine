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



void simulate_sand(Simulation* simulation, float dt) {

    Profiling::start_counter(0);

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
    std::vector<glm::vec3>& velocities = simulation->velocities;


    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;

    //integration
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] += simulation->gravity * dt;
        //        if (attract_flag) {
        //            glm::vec3 p_to_cp = character_pos - positions[i];
        ////        float force_magnitude = 10 * std::min(1.0, 1.0/ glm::dot(p_to_gs, p_to_gs));
        //            velocities[i] += glm::normalize(p_to_cp) * simulation->particleRadius * 100.0f * dt * w;
        //        }
        //        if (blow_flag) {
        //            glm::vec3 p_to_cp = character_pos - positions[i];
        ////        float force_magnitude = 10 * std::min(1.0, 1.0/ glm::dot(p_to_gs, p_to_gs));
        //            velocities[i] += -glm::normalize(p_to_cp) * simulation->particleRadius * 1.0f * w;
        //        }

        positions_star[i] = positions[i] + velocities[i] * dt; // update both
    }
    //glm::vec3 *positions_tmp = new glm::vec3[simulation->total_allocated];
    memcpy(simulation->positions_tmp, positions_star, sizeof(glm::vec3) * simulation->total_allocated);

    //find_neighbors_counting_sort(simulation);
    Profiling::start_counter(3);
    find_neighbors_uniform_grid(simulation);
    Profiling::stop_counter(3);
    //    find_neighbors_brute_force(simulation);

        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    for (int substep = 0; substep < 4; ++substep) {

        std::vector<std::pair<int, int>> contacts;
        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
            for (int j : neighbors[i]) {
                if (i >= j) continue;
                glm::vec3 ij = simulation->positions_tmp[i] - simulation->positions_tmp[j];
                float len = glm::length(ij);
                if (len < simulation->particleDiameter) {
                    contacts.push_back({ i, j });
                }
            }
        }

        // collision and friction
        for (auto& contact : contacts) {
            int i = contact.first, j = contact.second;
            glm::vec3 ij = simulation->positions_tmp[i] - simulation->positions_tmp[j];
            float len = glm::length(ij);

            // TODO: can be optimized by using a single mass;
            // i is always sand
            float wi = 1.0f / simulation->mass;
            // compute j's mass
            float wj = (j >= simulation->ptr_sand_start && j < simulation->ptr_sand_end) ? 1.0f / simulation->mass
                : 0.0f;
            float wi_wiPwj = wi / (wi + wj);
            float wj_wiPwj = wj / (wi + wj);

            glm::vec3 delta_ij = (len - simulation->particleDiameter) * ij / avoid0(len);
            glm::vec3 delta_pi = -wi_wiPwj * delta_ij;
            glm::vec3 delta_pj = +wj_wiPwj * delta_ij;
            glm::vec3 positions_star_i = simulation->positions_tmp[i] + collision_coeff * delta_pi;
            glm::vec3 positions_star_j = simulation->positions_tmp[j] + collision_coeff * delta_pj;


            float d = simulation->particleDiameter - len;
            glm::vec3 n = glm::normalize(positions_star_i - positions_star_j);
            // FIXME: Should use initial positions but currently very unstable
//            glm::vec3 delta_x_star_ij = (positions_star_i - positions[i]) - (positions_star_j - positions[j]);
            glm::vec3 delta_x_star_ij = (positions_star_i - simulation->positions_tmp[i]) - (positions_star_j - simulation->positions_tmp[j]);
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

            positions_star[i] += positions_star_i - simulation->positions_tmp[i];
            positions_star[j] += positions_star_j - simulation->positions_tmp[j];

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


    }



    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] = (positions_star[i] - positions[i]) / simulation->time_step;
        positions[i] = positions_star[i];
    }

    //TODO: this perturbation is for test, remove this in the future
    static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.01, 0.01, 0.01);
        perturbation = false;
    }

    Profiling::stop_counter(2);
    Profiling::stop_counter(0);
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

    dt = glm::clamp(dt, 0.001f, 0.01f);
    dt = (1.0f / 30.0f);
    simulation->time_step = dt;

    int n = simulation->num_particles;
    int X = simulation->domainX;
    int Y = simulation->domainY;
    int Z = simulation->domainZ;

    glm::vec3* positions = simulation->positions;
    glm::vec3* positions_star = simulation->positions_star;
    std::vector<glm::vec3>& velocities = simulation->velocities;


    std::vector<float>& lambdas = simulation->lambdas;
    std::vector<std::vector<int>>& neighbors = simulation->neighbors;

    //integration

    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] += simulation->gravity * dt;
        positions_star[i] = positions[i] + velocities[i] * dt; // update both
    }

    //find_neighbors_counting_sort(simulation);
    Profiling::start_counter(3);
    find_neighbors_uniform_grid(simulation);
    Profiling::stop_counter(3);

    float particleDiameter2 = simulation->particleDiameter * simulation->particleDiameter;
        // solve contact constraints(collision, friction), http://mmacklin.com/flex_eurographics_tutorial.pdf
    
    Profiling::start_counter(4);
    for (int substep = 0; substep < 4; ++substep) {

        std::vector<std::pair<int, int>> contacts;
        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
            for (int j : neighbors[i]) {
                if (i >= j) continue;
                glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[j];
                float len = glm::dot(ij, ij); //glm::length(ij);
                if (len < particleDiameter2) {
                    contacts.push_back({ i, j });
                }
            }
        }

        float mass_inv = 1.0f / simulation->mass;
        // collision and friction
        for (auto& contact : contacts) {
            int i = contact.first, j = contact.second;
            glm::vec3 ij = simulation->positions_star[i] - simulation->positions_star[j];
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
//            glm::vec3 delta_x_star_ij = (positions_star_i - positions[i]) - (positions_star_j - positions[j]);
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

        // particle-bounadry collision
        for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
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

        //memcpy(simulation->positions_tmp + simulation->ptr_sand_start, positions_star + simulation->ptr_sand_start, sizeof(glm::vec3) * simulation->num_sand_particles);


    }

    Profiling::stop_counter(4);


    float dt_inv = 1.0f / simulation->time_step;
    //update velocities and sync positions
    for (int i = simulation->ptr_sand_start; i < simulation->ptr_sand_end; i++) {
        velocities[i] = (positions_star[i] - positions[i]) * dt_inv;
        positions[i] = positions_star[i];
    }

    //TODO: this perturbation is for test, remove this in the future
    static bool perturbation = true;
    if (perturbation) {
        velocities[simulation->ptr_sand_start] += glm::vec3(0.01, 0.01, 0.01);
        perturbation = false;
    }

    Profiling::stop_counter(2);
}

}