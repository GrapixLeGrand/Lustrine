#pragma once

#include "BulletPhysicsSimulation.hpp"

namespace Lustrine {

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 */
	extern void init_bullet(BulletPhyicsSimulation* simulation);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 */
	extern void clean_bullet(BulletPhyicsSimulation* simulation);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param dt 
	 */
	extern void simulate_bullet(BulletPhyicsSimulation* simulation, float dt);

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param position 
	 * @param is_dynamic 
	 * @return int the index of the body in the main body array
	 */
	extern int add_box(BulletPhyicsSimulation* simulation, glm::vec3& position, bool is_dynamic, glm::vec4 color);

	//TODO
	extern void apply_force(BulletPhyicsSimulation* simulation, int body_index, glm::vec3 force);
	
	//TODO
	extern void set_velocity(BulletPhyicsSimulation* simulation, int body_index, glm::vec3 velocity);

}