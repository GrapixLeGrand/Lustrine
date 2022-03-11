
#include "BulletPhysics.hpp"
#include "BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h"

namespace Lustrine {

	void init_bullet(BulletPhyicsSimulation* simulation) {
		
		simulation->collisionConfiguration = new btDefaultCollisionConfiguration();
		simulation->dispatcher = new btCollisionDispatcher(simulation->collisionConfiguration);
		simulation->broadPhase = new btDbvtBroadphase();
		btGImpactCollisionAlgorithm::registerAlgorithm(simulation->dispatcher);

		//simulation->solver = new btSequentialImpulseConstraintSolver();
		//simulation->dynamicWorld = new btDiscreteDynamicsWorld(simulation->dispatcher, simulation->broadPhase, simulation->solver, simulation->collisionConfiguration);
		btVector3 gravity(0, -10, 0);
		simulation->dynamicWorld->setGravity(gravity);

	}
}