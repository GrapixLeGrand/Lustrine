
#include "BulletPhysics.hpp"
#include "BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h"
#include "glm/ext.hpp"
#include <iostream>


namespace Lustrine {
namespace Bullet {

	btVector3 glmToBullet(glm::vec3& v) {
		btVector3 result;
		result.setX(v.x);
		result.setY(v.y);
		result.setZ(v.z);
		return result;
	}

	glm::vec3 bulletToGlm(btVector3& v) {
		glm::vec3 result;
		result.x = v.getX();
		result.y = v.getY();
		result.z = v.getZ();
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
	 btRigidBody* bullet_create_rigidbody(Simulation* simulation, float mass, const btTransform& startTransform, btCollisionShape* shape, int index) {

		btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));
		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);
		btVector3 localInertia(0, 0, 0);
		btRigidBody* body = nullptr;
		if (isDynamic) {
			shape->calculateLocalInertia(mass, localInertia);
			//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);
			body = new btRigidBody(cInfo);
		
		} else {
			btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			body = new btRigidBody(mass, myMotionState, shape, localInertia);
			body->setActivationState(ISLAND_SLEEPING);
			body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_STATIC_OBJECT);
		}

		//body->setContactProcessingThreshold(m_defaultContactProcessingThreshold);
		
		//body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT); //warning was dynamic object before.
		//body->activate(true);

		//btRigidBody* body = new btRigidBody(mass, 0, shape, localInertia);
		//body->setWorldTransform(startTransform);
        //

		body->setUserIndex(index);
		return body;
	}

	void init_bullet(Simulation* simulation) {
		
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

		simulation->bodies_collisions = std::vector<std::vector<int>>({});
	}

	void clean_bullet(Simulation* simulation) {

		std::cout << "Lustrine::info Exiting bullet" << std::endl;

		delete simulation->dynamicWorld;
		delete simulation->solver;
		delete simulation->broadPhase;
		delete simulation->dispatcher;
		delete simulation->collisionConfiguration;
	}

	void simulate_bullet(Simulation* simulation, float dt) {
		gather_collisions(simulation);
		//print_collisions(simulation);
		simulation->dynamicWorld->stepSimulation(dt);
	}

	int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic) {
		return add_box(simulation, simulation->unit_box_shape, position, is_dynamic, btBroadphaseProxy::AllFilter, INT32_MAX);
	}

	int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, int group, int mask) {
		return add_box(simulation, simulation->unit_box_shape, position, is_dynamic, group, mask);
	}

	int add_box(Simulation* simulation, btBoxShape* box_shape, glm::vec3 position, bool is_dynamic, int group, int mask) {
		btTransform tmpTransform;
		tmpTransform.setIdentity();
		btVector3 pos = glmToBullet(position);
		tmpTransform.setOrigin(pos);
		float mass = is_dynamic ? simulation->box_mass : 0.0f;
		btRigidBody* newBody = bullet_create_rigidbody(simulation, mass, tmpTransform, box_shape, simulation->num_bodies);
		newBody->setFriction(simulation->default_body_friction);
		simulation->rigidbodies.push_back(newBody);
		simulation->transforms.push_back(tmpTransform);
		int result = simulation->num_bodies;
		simulation->num_bodies++;
		simulation->bodies_collisions.resize(simulation->num_bodies, {});
		simulation->dynamicWorld->addRigidBody(newBody, group, mask);
		return result;
	}

	int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims, int group, int mask) {
		
		btBoxShape* new_box_shape = new btBoxShape(glmToBullet(half_dims)); //for now we register only a unit cube
		simulation->collisionShapes.push_back(new_box_shape);
		
		int result = add_box(simulation, new_box_shape, position, is_dynamic, group, mask);
		return result;
		/*btTransform tmpTransform;
		tmpTransform.setIdentity();
		btVector3 pos = glmToBullet(position);
		tmpTransform.setOrigin(pos);
		float mass = is_dynamic ? simulation->box_mass : 0.0f;
		btRigidBody* newBody = bullet_create_rigidbody(simulation, mass, tmpTransform, new_box_shape, simulation->num_bodies);
		newBody->setFriction(1.0f);
		simulation->rigidbodies.push_back(newBody);
		simulation->transforms.push_back(tmpTransform);
		int result = simulation->num_bodies;
		simulation->num_bodies++;
		simulation->bodies_collisions.resize(simulation->num_bodies, {});
		simulation->dynamicWorld->addRigidBody(newBody, simulation->collision_group_0, simulation->collision_mask_0);
		return result;*/
	}

	void allocate_particles_colliders(Simulation* simulation, int num_particles) {
		
		if (num_particles < simulation->sand_particles_colliders.size()) {
			return;
		}

		size_t old_size = 0;
		if (simulation->sand_particles_colliders.size() != 0) {
			old_size = simulation->sand_particles_colliders.size();
			simulation->sand_particles_colliders.resize(num_particles - old_size, nullptr);
		} else {
			simulation->sand_particles_colliders = std::vector<btRigidBody*>(num_particles, nullptr);
		}

		simulation->bodies_collisions.resize(simulation->num_bodies + num_particles, std::vector<int>{});

		for (old_size; old_size < num_particles; old_size++) {
			btTransform tmpTransform;
			tmpTransform.setIdentity();
			btVector3 pos (0.0f, 0.0f, 0.0f);
			tmpTransform.setOrigin(pos);
			float mass = 0.0f; //particles must not be moving at first
			simulation->sand_particles_colliders[old_size] = bullet_create_rigidbody(simulation, mass, tmpTransform, simulation->unit_box_shape, simulation->num_bodies);
			simulation->dynamicWorld->addRigidBody(simulation->sand_particles_colliders[old_size], simulation->collision_group_1, simulation->collision_mask_1);
			simulation->num_bodies++;
		}
	}

	void set_body_no_rotation(Simulation* simulation, int body_index) {
		btVector3 linFact (0.0, 0.0, 0.0);
        simulation->rigidbodies[body_index]->setAngularFactor(linFact);
	}

	void set_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity) {
		simulation->rigidbodies[body_index]->activate(true);
        simulation->rigidbodies[body_index]->setLinearVelocity(glmToBullet(velocity));
	}

	void print_resume(const Simulation* simulation) {
		std::cout << "Lustrine::Bullet " << "\n" 
			<< "\tnum registered bodies: " << simulation->num_bodies << "\n"
			<< "\tnum registered shapes: " << simulation->collisionShapes.size() << std::endl;
	}

	void apply_impulse(Simulation* simulation, int body_index, glm::vec3 impulse, glm::vec3 position) {
		simulation->rigidbodies[body_index]->applyImpulse(glmToBullet(impulse), glmToBullet(position));
		
	}

	void print_collisions(Simulation* simulation) {
		for (int i = 0; i < simulation->bodies_collisions.size(); i++) {
			std::cout << i << " -> " << std::endl;
			for (int j = 0; j < simulation->bodies_collisions[i].size(); j++) {
				std::cout << simulation->bodies_collisions[i][j] << ", ";
			}
			std::cout << std::endl;
		}
	}

	bool check_collision(Simulation* simulation, int body1, int body2) {
		for (int i = 0; i < simulation->bodies_collisions[body1].size(); i++) {
			if (simulation->bodies_collisions[body1][i] == body2) {
				return true;
			}
		}
		return false;
	}

	void gather_collisions(Simulation* simulation) {

		for (int i = 0; i < simulation->bodies_collisions.size(); i++) {
			simulation->bodies_collisions[i].clear();
		}

		btDispatcher* dp = simulation->dispatcher;// world->getDispatcher();
    	const int numManifolds = dp->getNumManifolds();
		for (int m=0; m<numManifolds; ++m) {
            
			btPersistentManifold* man = dp->getManifoldByIndexInternal(m);
            const btRigidBody* obA = static_cast<const btRigidBody*>(man->getBody0());
            const btRigidBody* obB = static_cast<const btRigidBody*>(man->getBody1());

			int indexA = obA->getUserIndex();
			int indexB = obB->getUserIndex();

			simulation->bodies_collisions[indexA].push_back(indexB);
			simulation->bodies_collisions[indexB].push_back(indexA);	

    	}

	}

	void set_particles_box_colliders_positions(Simulation* simulation, glm::vec3* particles, int start, int end) {
		
		for (int i = 0; i < simulation->sand_particles_colliders.size(); i++) {
			btTransform& t = simulation->sand_particles_colliders[i]->getWorldTransform();
			t.setOrigin(glmToBullet(particles[i]));
        	simulation->sand_particles_colliders[i]->getMotionState()->setWorldTransform(t);
		}

	}

	void set_particles_positions(Simulation* simulation, int body, std::vector<glm::vec3> particles, int start, int end, float particleRadius) {
		
		btTransform& t = simulation->rigidbodies[body]->getWorldTransform();
		glm::mat4 m (0.0f);
		t.getOpenGLMatrix(glm::value_ptr(m));
		glm::vec3 box_dims (0.0); 
		box_dims *= 2;

		float particleDiameter = 2.0f * particleRadius;
		int X = (int) (box_dims.x / particleDiameter);
		int Y = (int) (box_dims.y / particleDiameter);
		int Z = (int) (box_dims.z / particleDiameter);

		int num_particles_to_setup = end - start;
		assert(X * Y * Z == num_particles_to_setup);

		glm::vec3 offset = bulletToGlm(t.getOrigin());
		
		for (int x = 0; x < X; x++) {
			for (int y = 0; y < Y; y++) {
				for (int z = 0; z < Z; z++) {
					
				}
			}
		}


	}
}
}