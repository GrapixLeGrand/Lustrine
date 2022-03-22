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
		btVector3 gravity;

		btAlignedObjectArray<btCollisionShape*> collisionShapes;
		btBoxShape* unit_box_shape; //unit box shape
		
		int num_bodies = 0;

		std::vector<btTransform> transforms;
		std::vector<btRigidBody*> rigidbodies; //store all registered bodies
		std::vector<std::vector<int>> bodies_collisions;

		std::vector<btRigidBody*> sand_particles_colliders;

		//attention this is work in progress and this group thing does not work for now
		int collision_group_0 = (1 << 11);
		int collision_mask_0 = (1 << 11); //player + objects
		int collision_group_1 = (1 << 12); //particles (and optionally player)
		int collision_mask_1 = (1 << 12);
		int collision_group_2 = (1 << 13); //particles (and optionally player)
		int collision_mask_2 = (1 << 13);

		float default_body_friction = 5.0f;
		float box_mass = 1.0f;

	};

}
}