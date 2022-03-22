#pragma once

#include "BulletPhysicsSimulation.hpp"

namespace Lustrine {
namespace Bullet {

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 */
	extern void init_bullet(Simulation* simulation);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 */
	extern void clean_bullet(Simulation* simulation);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param dt 
	 */
	extern void simulate_bullet(Simulation* simulation, float dt);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param position 
	 * @param is_dynamic 
	 * @return int 
	 */
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic);
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, int group, int mask);
	extern int add_box(Simulation* simulation, btBoxShape* box_shape, glm::vec3 position, bool is_dynamic, int group, int mask);
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, int group, int mask);
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims, int group, int mask);

	/**
	 * @brief Allocates a pool of bodies that could collide with the rest of the world
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern void allocate_particles_colliders(Simulation* simulation, int num_particles);

	/**
	 * @brief Body will move but not rotate (main player)
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern void set_body_no_rotation(Simulation* simulation, int body_index);
	
	/**
	 * @brief Set the velocity of the body
	 * 
	 * @param simulation 
	 * @param body_index 
	 * @param velocity 
	 */
	extern void set_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity);

	/**
	 * @brief Apply an impulse on the center of mass of the body
	 * 
	 * @param simulation 
	 * @param body_index 
	 * @param force 
	 */
	extern void apply_impulse(Simulation* simulation, int body_index, glm::vec3 impulse, glm::vec3 relative_position);


	/**
	 * @brief Clear and fill the collisions for all the bodies that collided (this function is called in the simulate bullet)
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern void gather_collisions(Simulation* simulation);

	/**
	 * @brief returns true if the body considered do collides (are touching in the current frame), false otw.
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern bool check_collision(Simulation* simulation, int body1, int body2);


	/**
	 * @brief helper to print the collisions
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern void print_collisions(Simulation* simulation);

	extern void set_particles_positions(Simulation* simulation, int body, std::vector<glm::vec3> particles, int start, int end, float particleRadius);

	//TODO
	extern void apply_force(Simulation* simulation, int body_index, glm::vec3 force);
	

	extern void print_resume(const Simulation* simulation);
}
}