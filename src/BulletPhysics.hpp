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
	//extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic); //, int group, int mask);
	extern int add_shape(Simulation* simulation, btConvexShape* box_shape, glm::vec3 position, bool is_dynamic); //, int group, int mask);
	//extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, int group, int mask);
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims);//, int group, int mask);


	extern int add_capsule(Simulation* simulation, glm::vec3 position);

	/**
	 * @brief set the gravity of bullet world
	 * 
	 * @param simulation 
	 * @param new_gravity 
	 */
	extern void set_gravity(Simulation* simulation, glm::vec3 new_gravity);
	extern glm::vec3 get_gravity(Simulation* simulation);

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
	 * @brief store in indices the indices of body colliding with body. Warning, the
	 * pointer is assume to have at least num_bodies entries.
	 * 
	 * @param simulation 
	 * @param body 
	 * @param indices 
	 * @param size 
	 */
	extern void check_collisions(Simulation* simulation, int body, int* indices, int* size);

	/**
	 * @brief Get the num bodies object
	 * 
	 * @param simulation 
	 * @return int 
	 */
	extern int get_num_bodies(Simulation* simulation);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param body 
	 * @return true 
	 * @return false 
	 */
	extern bool do_collide(Simulation* simulation, int body);

	/**
	 * @brief Set the particles box colliders positions object
	 * 
	 * @param simulation 
	 * @param particles 
	 * @param start 
	 * @param end 
	 */
	extern void set_particles_box_colliders_positions(Simulation* simulation, glm::vec3* particles, int start, int end);
	
	/**
	 * @brief returns the current body positions as specified by its motion state.
	 * 
	 * @param simulation 
	 * @param body 
	 * @return glm::vec3 
	 */
	extern glm::vec3 get_body_position(Simulation* simulation, int body);

	/**
	 * @brief Will override the motion state of the object and set its position to be the new
	 * one specified
	 * 
	 * @param simulation 
	 * @param body 
	 * @param new_position 
	 */
	extern void set_body_position(Simulation* simulation, int body, glm::vec3 new_position);

	/**
	 * @brief helper to print the collisions
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	extern void print_collisions(Simulation* simulation);
	extern void set_particles_positions(Simulation* simulation, int body, std::vector<glm::vec3> particles, int start, int end, float particleRadius);
	extern void apply_force(Simulation* simulation, int body_index, glm::vec3 force);
	

	extern void print_resume(const Simulation* simulation);
}
}