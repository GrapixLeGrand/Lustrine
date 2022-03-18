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

	struct BulletPhyicsSimulation {

		btBroadphaseInterface* broadPhase;
		btCollisionDispatcher* dispatcher;
		btConstraintSolver* solver;
		btDefaultCollisionConfiguration* collisionConfiguration;
		btDiscreteDynamicsWorld* dynamicWorld;
		btVector3 gravity;

		btAlignedObjectArray<btCollisionShape*> collisionShapes;
		btBoxShape* unit_box_shape; //unit box shape
		
		int num_bodies = 0;
		//int ptr_dynamic_start = 0;
		//int ptr_dynamic_end = 0;
		//int ptr_static_start = 0;
		//int ptr_static_end = 0;

		std::vector<btTransform> transforms;
		std::vector<btRigidBody*> rigidbodies; //store all registered bodies
		//std::vector<glm::vec3> positions; //stores the positions of ALL boxes
		//std::vector<glm::vec4> colors;

		float box_mass = 1.0f;

	};

}