
#include "BulletPhysics.hpp"
#include "BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h"
#include <iostream>


namespace Lustrine {

	btVector3 glmToBullet(glm::vec3& v) {
		btVector3 result;
		result.setX(v.x);
		result.setY(v.y);
		result.setZ(v.z);
		return result;
	}

	/**
	 * @brief Create a Rigid Body object
	 * helper method taken from past project to generate a rigid body.
	 * The rigidbody is added to the simulation even if returned.
	 * @param mass if zero then body is static
	 * @param startTransform 
	 * @param shape 
	 * @return btRigidBody* 
	 */
	 btRigidBody* bullet_create_rigidbody(BulletPhyicsSimulation* simulation, float mass, const btTransform& startTransform, btCollisionShape* shape) {

		btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));
		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);
		btVector3 localInertia(0, 0, 0);
		if (isDynamic)
			shape->calculateLocalInertia(mass, localInertia);
			//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects

        #define USE_MOTIONSTATE 1
        #ifdef USE_MOTIONSTATE
		btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);

		btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);

		btRigidBody* body = new btRigidBody(cInfo);
		//body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);
		
		body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT); //warning was dynamic object before.
		//body->activate(true);

        #else
		btRigidBody* body = new btRigidBody(mass, 0, shape, localInertia);
		body->setWorldTransform(startTransform);
        #endif  //

		body->setUserIndex(-1);
		simulation->dynamicWorld->addRigidBody(body);
		return body;
	}

	void init_bullet(BulletPhyicsSimulation* simulation) {
		
		std::cout << "Lustrine::info Initializing bullet" << std::endl;

		simulation->collisionConfiguration = new btDefaultCollisionConfiguration();
		simulation->dispatcher = new btCollisionDispatcher(simulation->collisionConfiguration);
		simulation->broadPhase = new btDbvtBroadphase();
		btGImpactCollisionAlgorithm::registerAlgorithm(simulation->dispatcher);

		simulation->solver = new btSequentialImpulseConstraintSolver();
		simulation->dynamicWorld = new btDiscreteDynamicsWorld(simulation->dispatcher, simulation->broadPhase, simulation->solver, simulation->collisionConfiguration);
		btVector3 gravity(0, -10, 0);
		simulation->dynamicWorld->setGravity(gravity);

		btVector3 half(0.5, 0.5, 0.5);
		simulation->unit_box_shape = new btBoxShape(half); //for now we register only a unit cube
		simulation->collisionShapes.push_back(simulation->unit_box_shape);

	}

	void clean_bullet(BulletPhyicsSimulation* simulation) {

		std::cout << "Lustrine::info Exiting bullet" << std::endl;

		delete simulation->dynamicWorld;
		delete simulation->solver;
		delete simulation->broadPhase;
		delete simulation->dispatcher;
		delete simulation->collisionConfiguration;
	}

	void simulate_bullet(BulletPhyicsSimulation* simulation, float dt) {
		simulation->dynamicWorld->stepSimulation(dt);
	}

	int add_box(BulletPhyicsSimulation* simulation, glm::vec3& position, bool is_dynamic, glm::vec4 color = glm::vec4 {1.0, 0.0, 0.0, 1.0}) {
		btTransform tmpTransform;
		tmpTransform.setIdentity();
		btVector3 pos = glmToBullet(position);
		tmpTransform.setOrigin(pos);
		float mass = is_dynamic ? simulation->box_mass : 0.0f;
		btRigidBody* newBody = bullet_create_rigidbody(simulation, mass, tmpTransform, simulation->unit_box_shape);
		simulation->rigidbodies.push_back(newBody);
		simulation->positions.push_back(position);
		simulation->transforms.push_back(tmpTransform);
		int result = simulation->num_bodies;
		simulation->num_bodies++;
		simulation->colors.push_back(color);
		return result;
	}

}