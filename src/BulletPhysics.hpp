#pragma once

#include "Lustrine_Export.h"
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
	
	extern LUSTRINE_EXPORT void init_bullet(Simulation* simulation);
	extern LUSTRINE_EXPORT void clean_bullet(Simulation* simulation);
	extern LUSTRINE_EXPORT void simulate_bullet(Simulation* simulation, float dt, int start, int end);
	extern LUSTRINE_EXPORT int get_num_bodies(Simulation* simulation);
	extern LUSTRINE_EXPORT void set_gravity(Simulation* simulation, glm::vec3 new_gravity);
	extern LUSTRINE_EXPORT glm::vec3 get_gravity(Simulation* simulation);

	/**
	 * @brief Shapes creations
	 */

	extern LUSTRINE_EXPORT int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic);
	extern LUSTRINE_EXPORT int add_shape(Simulation* simulation, btConvexShape* box_shape, glm::vec3 position, bool is_dynamic); //, int group, int mask);
	extern LUSTRINE_EXPORT int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims);//, int group, int mask);
	
	extern LUSTRINE_EXPORT int add_detector_block(Simulation* simulation, glm::vec3 position, glm::vec3 half_dims);
	extern LUSTRINE_EXPORT int add_detector_cylinder(Simulation* simulation, glm::vec3 position, glm::vec3 half_dims);
	extern LUSTRINE_EXPORT int add_capsule(Simulation* simulation, glm::vec3 position, float radius, float height);

	
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

	extern LUSTRINE_EXPORT void apply_impulse(Simulation* simulation, int body_index, glm::vec3 impulse, glm::vec3 relative_position);
	//extern void apply_force(Simulation* simulation, int body_index, glm::vec3 force);
	extern LUSTRINE_EXPORT glm::vec3 get_body_position(Simulation* simulation, int body);
	extern LUSTRINE_EXPORT void set_body_position(Simulation* simulation, int body, glm::vec3 position);
	extern LUSTRINE_EXPORT void set_body_position(Simulation* simulation, int body, glm::vec3 new_position);
	extern LUSTRINE_EXPORT void set_body_rotations(Simulation* simulation, int body_index, bool X, bool Y, bool Z);
	extern LUSTRINE_EXPORT void set_body_no_rotation(Simulation* simulation, int body_index);
	extern LUSTRINE_EXPORT glm::vec3 get_body_velocity(Simulation* simulation, int body_index);
	extern LUSTRINE_EXPORT void set_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity);
	extern LUSTRINE_EXPORT void add_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity);
	extern LUSTRINE_EXPORT void set_body_frixion(Simulation* simulation, int body, float frixion);
	extern LUSTRINE_EXPORT float get_body_frixion(Simulation* simulation, int body);
	extern LUSTRINE_EXPORT void set_body_damping(Simulation* simulation, int body, float linear, float angular);
	extern LUSTRINE_EXPORT float get_body_lin_damping(Simulation* simulation, int body);
	//extern LUSTRINE_EXPORT float get_body_angular_damping(Simulation* simulation, int body)

	/**
	 * @brief Collision related 
	 */
	
	extern LUSTRINE_EXPORT void gather_collisions(Simulation* simulation);
	extern LUSTRINE_EXPORT bool check_collision(Simulation* simulation, int body1, int body2);
	extern LUSTRINE_EXPORT void check_collisions(Simulation* simulation, int body, int* indices, int* size);
	extern LUSTRINE_EXPORT bool do_collide(Simulation* simulation, int body);
	extern LUSTRINE_EXPORT bool do_collide_except_for(Simulation* simulation, int body, int exception_id);

	/**
	 * @brief Particles/rigidbodies related function
	 */

	extern LUSTRINE_EXPORT void allocate_particles_colliders(Simulation* simulation, int num_particles, float radius);
	extern LUSTRINE_EXPORT void set_particles_box_colliders_positions(Simulation* simulation, glm::vec3* particles, int start_ptr, int end_ptr);
	extern LUSTRINE_EXPORT void disable_particles_bounding_boxes(Simulation* simulation);
	extern LUSTRINE_EXPORT void enable_particles_bounding_boxes(Simulation* simulation);
	extern LUSTRINE_EXPORT void bind_foreign_sand_positions(Simulation* simulation, glm::vec3* foreign_positions);


	/**
	 * @brief Helpers
	 */

	extern LUSTRINE_EXPORT void print_resume(const Simulation* simulation);
	extern LUSTRINE_EXPORT void print_collisions(Simulation* simulation);
}
}