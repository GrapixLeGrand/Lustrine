#pragma once


#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
#endif

#include "btBulletDynamicsCommon.h"
#include <vector>

namespace Lustrine {
namespace Bullet {

	struct Simulation {

		btBroadphaseInterface* broadPhase;
		btCollisionDispatcher* dispatcher;
		btConstraintSolver* solver;
		btDefaultCollisionConfiguration* collisionConfiguration;
		btDiscreteDynamicsWorld* dynamicWorld;
		btVector3 gravity = btVector3(0.0f, -25.0f, 0.0f);

		btAlignedObjectArray<btCollisionShape*> collisionShapes;//seems to be auto managed memory
		btBoxShape* unit_box_shape; //unit box shape
		
		int num_bodies = 0;
		int num_particles_allocated = 100;//number of boxes colliders around the player
		int player_id = 0;
		glm::vec3 player_position = glm::vec3(0.0f);
		float player_box_radius = 4.0f;//radius up to which we try to add boxes


		bool allocated_particles_bounding_boxes = false;//FOR NOW A SINGLE ALLOCATION IS POSSIBLE
		bool particles_bounding_box_requested_state = true;//By default the particles are enabled
		bool particles_bounding_box_current_state = true;//internal variable used to detect rising edge on bounding box deactivation/activation

		glm::vec3* foreign_sand_positions = nullptr; //ptr not owned
		size_t ptr_bounding_box_start = 0;//inclusive
		size_t ptr_bounding_box_end = 0;//exclusive

		std::vector<btTransform> transforms;
		std::vector<btRigidBody*> rigidbodies; //store all registered bodies
		std::vector<std::vector<int>> bodies_collisions;


		//std::vector<btRigidBody*> sand_particles_colliders;

		//attention this is work in progress and this group thing does not work for now
		int collision_group_0 = (1 << 11);
		int collision_mask_0 = (1 << 11); //player + objects
		int collision_group_1 = (1 << 12); //particles (and optionally player)
		int collision_mask_1 = (1 << 12);
		int collision_group_2 = (1 << 13); //particles (and optionally player)
		int collision_mask_2 = (1 << 13);

		float default_body_friction = 0.0f;
		float box_mass = 1.0f;

	};

}
}