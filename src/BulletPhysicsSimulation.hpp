#pragma once


#ifdef PLATFORM_WINDOWS
#include "../thirdparty/glm-0.9.9.8/glm/glm.hpp"
#endif
#ifdef PLATFORM_UNIX
#include "glm/glm.hpp"
#endif

#include "btBulletDynamicsCommon.h"

namespace Lustrine {

	struct BulletPhyicsSimulation {

		btCollisionShape* collisionShapes;
		btBroadphaseInterface* broadPhase;
		btCollisionDispatcher* dispatcher;
		btConstraintSolver* solver;
		btDefaultCollisionConfiguration* collisionConfiguration;
		btDiscreteDynamicsWorld* dynamicWorld;
		btVector3 gravity;

	};

}