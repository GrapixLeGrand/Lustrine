#pragma once

#include "BulletPhysicsSimulation.hpp"

namespace Lustrine {
namespace Bullet {

	/**
	 * @brief physical types of the bodies
	 * 
	 */
	enum BodyType {
		DYNAMIC,//body that moves
		KINEMATIC,//not moving but movable by user.
		DETECTOR//not moving and not interacting but collision happens.
	};
	
	extern void init_bullet(Simulation* simulation);
	extern void clean_bullet(Simulation* simulation);
	extern void simulate_bullet(Simulation* simulation, float dt);
	extern int get_num_bodies(Simulation* simulation);
	extern void set_gravity(Simulation* simulation, glm::vec3 new_gravity);
	extern glm::vec3 get_gravity(Simulation* simulation);

	/**
	 * @brief Shapes creations
	 */

	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic);
	extern int add_shape(Simulation* simulation, btConvexShape* box_shape, glm::vec3 position, bool is_dynamic); //, int group, int mask);
	extern int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims);//, int group, int mask);
	
	extern int add_detector_block(Simulation* simulation, glm::vec3 position, glm::vec3 half_dims);
	extern int add_capsule(Simulation* simulation, glm::vec3 position, float radius, float height);

	
	/**
	 * @brief clear the velocity of the body
	 * TODO impl
	 * @param simulation 
	 * @param body_index 
	 * @param velocity 
	 */
	//extern void clear_body_velocity(Simulation* simulation, int body_index);

	/**
	 * @brief bodies properties setters
	 */

	extern void apply_impulse(Simulation* simulation, int body_index, glm::vec3 impulse, glm::vec3 relative_position);
	extern void apply_force(Simulation* simulation, int body_index, glm::vec3 force);
	extern glm::vec3 get_body_position(Simulation* simulation, int body);
	extern void set_body_position(Simulation* simulation, int body, glm::vec3 new_position);
	extern void set_body_rotations(Simulation* simulation, int body_index, bool X, bool Y, bool Z);
	extern void set_body_no_rotation(Simulation* simulation, int body_index);
	extern void set_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity);
	extern void add_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity);

	/**
	 * @brief Collision related 
	 */
	
	extern void gather_collisions(Simulation* simulation);
	extern bool check_collision(Simulation* simulation, int body1, int body2);
	extern void check_collisions(Simulation* simulation, int body, int* indices, int* size);
	extern bool do_collide(Simulation* simulation, int body);

	
	
	
	extern void allocate_particles_colliders(Simulation* simulation, int num_particles);
	extern void set_particles_box_colliders_positions(Simulation* simulation, glm::vec3* particles, int start, int end);
	extern void set_particles_positions(Simulation* simulation, int body, std::vector<glm::vec3> particles, int start, int end, float particleRadius);
	extern void disable_particles_bounding_boxes(Simulation* simulation);
		
	/**
	 * @brief Helpers
	 */

	extern void print_resume(const Simulation* simulation);
	extern void print_collisions(Simulation* simulation);
}
}