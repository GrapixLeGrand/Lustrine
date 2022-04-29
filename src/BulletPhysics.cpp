
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

	glm::vec3 bulletToGlm(const btVector3& v) {
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
	 btRigidBody* bullet_create_rigidbody(Simulation* simulation, BodyType type, float mass, const btTransform& startTransform, btCollisionShape* shape, int index) {

		btAssert((!shape || shape->getShapeType() != INVALID_SHAPE_PROXYTYPE));
		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);
		btVector3 localInertia(0, 0, 0);
		btRigidBody* body = nullptr;
		btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
		if (type == DYNAMIC) {
			shape->calculateLocalInertia(mass, localInertia);
			//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
			btRigidBody::btRigidBodyConstructionInfo cInfo(mass, myMotionState, shape, localInertia);
			body = new btRigidBody(cInfo);
		} else if (type == KINEMATIC) {
			//btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			body = new btRigidBody(mass, myMotionState, shape, localInertia);
			body->setActivationState(ISLAND_SLEEPING);
			body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_KINEMATIC_OBJECT);// CF_STATIC_OBJECT);
		} else if (type == DETECTOR) {
			//btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
			body = new btRigidBody(mass, myMotionState, shape, localInertia);
			body->setActivationState(ISLAND_SLEEPING);
			body->setCollisionFlags(body->getCollisionFlags() | btCollisionObject::CF_STATIC_OBJECT | btCollisionObject::CF_NO_CONTACT_RESPONSE);
		} else {	
			assert(false);
		}

		body->setUserIndex(index);
		return body;
	}

	/**
	 * @brief initialize the engine. Must be called before anything regarding to bullet.
	 * 
	 * @param simulation 
	 */
	void init_bullet(Simulation* simulation) {
		
		std::cout << "Lustrine::info Initializing bullet" << std::endl;

		simulation->collisionConfiguration = new btDefaultCollisionConfiguration();
		simulation->dispatcher = new btCollisionDispatcher(simulation->collisionConfiguration);
		simulation->broadPhase = new btDbvtBroadphase();
		btGImpactCollisionAlgorithm::registerAlgorithm(simulation->dispatcher);

		simulation->solver = new btSequentialImpulseConstraintSolver();
		simulation->dynamicWorld = new btDiscreteDynamicsWorld(simulation->dispatcher, simulation->broadPhase, simulation->solver, simulation->collisionConfiguration);
		//btVector3 gravity(0, -10, 0);
		simulation->dynamicWorld->setGravity(simulation->gravity);

		btVector3 half(0.5, 0.5, 0.5);
		simulation->unit_box_shape = new btBoxShape(half); //for now we register only a unit cube
		simulation->collisionShapes.push_back(simulation->unit_box_shape);

		simulation->bodies_collisions = std::vector<std::vector<int>>({});
	}

	/**
	 * @brief set the gravity of bullet world
	 * 
	 * @param simulation 
	 * @param new_gravity 
	 */
	void set_gravity(Simulation* simulation, glm::vec3 new_gravity) {
		simulation->dynamicWorld->setGravity(glmToBullet(new_gravity));
	}

	/**
	 * @brief get the world's gravity
	 * 
	 * @param simulation 
	 * @return glm::vec3 
	 */
	glm::vec3 get_gravity(Simulation* simulation) {
		btVector3 gravity = simulation->dynamicWorld->getGravity();
		return bulletToGlm(gravity);
	}

	/**
	 * @brief destroy the memory allocated by the engine.
	 * 
	 * @param simulation 
	 */
	void clean_bullet(Simulation* simulation) {

		std::cout << "Lustrine::info Exiting bullet" << std::endl;

		delete simulation->dynamicWorld;
		delete simulation->solver;
		delete simulation->broadPhase;
		delete simulation->dispatcher;
		delete simulation->collisionConfiguration;

		for (int i = 0; i < simulation->rigidbodies.size(); i++) {
			delete simulation->rigidbodies[i]->getMotionState();
			delete simulation->rigidbodies[i];
		}
		simulation->rigidbodies.clear();
		simulation->num_bodies = 0;

		for (int i = 0; i < simulation->collisionShapes.size(); i++) {
			delete simulation->collisionShapes[i];
		}
		simulation->collisionShapes.clear();

		simulation->allocated_particles_bounding_boxes = false; // in case we reuse the simulation after destroyed
	
	}

	/**
	 * @brief 
	 * 
	 * @param simulation 
	 * @param dt 
	 */
	void simulate_bullet(Simulation* simulation, float dt, int sand_start, int sand_end) {
		gather_collisions(simulation);
		if (simulation->particles_bounding_box_current_state == false && simulation->particles_bounding_box_requested_state == true) {
			Lustrine::Bullet::enable_particles_bounding_boxes(simulation);
			simulation->particles_bounding_box_current_state = true;
		} else if (simulation->particles_bounding_box_current_state == true && simulation->particles_bounding_box_requested_state == false) {
			Lustrine::Bullet::disable_particles_bounding_boxes(simulation);
			simulation->particles_bounding_box_current_state = false;
		}
		set_particles_box_colliders_positions(simulation, simulation->foreign_sand_positions, sand_start, sand_end);
		simulation->dynamicWorld->stepSimulation(dt);
	}

	/**
	 * @brief adds a capsule at the given location. The capsule is necessrily dynamic and will fall.
	 * Prefered collider for player for now.
	 * 
	 * @param simulation 
	 * @param position 
	 * @param radius 
	 * @param height 
	 * @return int 
	 */
	int add_capsule(Simulation* simulation, glm::vec3 position, float radius, float height) {
		btCapsuleShape* capsule = new btCapsuleShape(radius, height);
		simulation->collisionShapes.push_back(capsule);
		int result = add_shape(simulation, capsule, position, true); //, group, mask);
		return result;
	}

	/**
	 * @brief adds a box and returns its body index
	 * 
	 * @param simulation 
	 * @param position 
	 * @param is_dynamic 
	 * @return int 
	 */
	int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic) {
		return add_shape(simulation, simulation->unit_box_shape, position, is_dynamic); //, btBroadphaseProxy::AllFilter, INT32_MAX);
	}

	int add_shape(Simulation* simulation, btConvexShape* box_shape, glm::vec3 position, bool is_dynamic) { //, int group, int mask) {
		btTransform tmpTransform;
		tmpTransform.setIdentity();
		btVector3 pos = glmToBullet(position);
		tmpTransform.setOrigin(pos);
		float mass = is_dynamic ? simulation->box_mass : 0.0f;
		BodyType type = mass > 0.0 ? DYNAMIC : KINEMATIC;
		btRigidBody* newBody = bullet_create_rigidbody(simulation, type, mass, tmpTransform, box_shape, simulation->num_bodies);
		newBody->setFriction(simulation->default_body_friction);
		//newBody->setLinearFactor(btVector3(1.5f, 1.5f, 1.5f));
		//newBody->setDamping(-100.0f, 0.0f);
		
		//btVector3 linearFactor = btVector3(0, 1.f, 0);
		//newBody->setLinearFactor(linearFactor);

		simulation->rigidbodies.emplace_back(newBody);
		simulation->transforms.push_back(tmpTransform);
		int result = simulation->num_bodies;
		simulation->num_bodies++;
		simulation->bodies_collisions.resize(simulation->num_bodies, {});
		simulation->dynamicWorld->addRigidBody(newBody);//, group, mask);
		return result;
	}

	int add_box(Simulation* simulation, glm::vec3 position, bool is_dynamic, glm::vec3 half_dims) { //, int group, int mask) {
		
		btBoxShape* new_box_shape = new btBoxShape(glmToBullet(half_dims)); //for now we register only a unit cube
		simulation->collisionShapes.push_back(new_box_shape);
		
		int result = add_shape(simulation, new_box_shape, position, is_dynamic); //, group, mask);
		return result;
	}

	/**
	 * @brief Add a ghost object with box shape. No collision response is generated, however collisions
	 * get registered.
	 * 
	 * @param simulation 
	 * @param position 
	 * @param is_dynamic 
	 * @param half_dims 
	 * @return int 
	 */
	int add_detector_block(Simulation* simulation, glm::vec3 position, glm::vec3 half_dims) {
		
		btTransform tmpTransform;
		tmpTransform.setIdentity();
		btVector3 pos = glmToBullet(position);
		tmpTransform.setOrigin(pos);

		btBoxShape* shape = new btBoxShape(glmToBullet(half_dims));
		simulation->collisionShapes.push_back(shape);

		btRigidBody* newBody = bullet_create_rigidbody(simulation, DETECTOR, 0.0f, tmpTransform, shape, simulation->num_bodies);
		newBody->setFriction(simulation->default_body_friction);
		
		simulation->rigidbodies.push_back(newBody);
		simulation->transforms.push_back(tmpTransform);
		int result = simulation->num_bodies;
		simulation->num_bodies++;
		simulation->bodies_collisions.resize(simulation->num_bodies, {});
		simulation->dynamicWorld->addRigidBody(newBody);//, group, mask);
		return result;
	}


	/**
	 * @brief Allocates a pool of bodies that could collide with the rest of the world
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
<<<<<<< HEAD
	void allocate_particles_colliders(Simulation* simulation, int num_particles, float radius) {
=======
	void allocate_particles_colliders(Simulation* simulation) {
>>>>>>> cec09dcf7d0348d9586ef58870ab02c8699b8aad
		
		assert(simulation->allocated_particles_bounding_boxes == false);
		if (simulation->allocated_particles_bounding_boxes == true) {
			std::cout << "already allocated particles! (not allowed to allocate twice" << std::endl;
		}

		simulation->bodies_collisions.resize(simulation->num_bodies + simulation->num_particles_allocated, std::vector<int>{});
		simulation->rigidbodies.resize(simulation->num_bodies + simulation->num_particles_allocated, nullptr);

		simulation->ptr_bounding_box_start = simulation->num_bodies;

		for (size_t old_size = 0; old_size < simulation->num_particles_allocated; old_size++) {
			btTransform tmpTransform;
			tmpTransform.setIdentity();
			btVector3 pos (0.0f, 0.0f, 0.0f);
			tmpTransform.setOrigin(pos);
			float mass = 0.0f; //particles must not be moving at first
			size_t index = old_size + simulation->num_bodies;

			simulation->rigidbodies[index] = bullet_create_rigidbody(simulation, KINEMATIC, mass, tmpTransform, new btSphereShape(radius), index);
			//simulation->rigidbodies[index] = bullet_create_rigidbody(simulation, KINEMATIC, mass, tmpTransform, new btBoxShape(btVector3(radius, radius, radius)), index);
			simulation->dynamicWorld->addRigidBody(simulation->rigidbodies[index]); //, simulation->collision_group_1, simulation->collision_mask_1);
		}

		simulation->num_bodies += simulation->num_particles_allocated;
		simulation->ptr_bounding_box_end = simulation->num_bodies;
		simulation->allocated_particles_bounding_boxes = true;

		assert(simulation->ptr_bounding_box_end - simulation->ptr_bounding_box_start == simulation->num_particles_allocated);

		//simulation->sand_particles_colliders = std::vector<btRigidBody*>(num_particles, nullptr);

		/*
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

		for (; old_size < num_particles; old_size++) {
			btTransform tmpTransform;
			tmpTransform.setIdentity();
			btVector3 pos (0.0f, 0.0f, 0.0f);
			tmpTransform.setOrigin(pos);
			float mass = 0.0f; //particles must not be moving at first
			simulation->sand_particles_colliders[old_size] = bullet_create_rigidbody(simulation, KINEMATIC, mass, tmpTransform, simulation->unit_box_shape, simulation->num_bodies);
			simulation->dynamicWorld->addRigidBody(simulation->sand_particles_colliders[old_size]); //, simulation->collision_group_1, simulation->collision_mask_1);
			simulation->num_bodies++;
		}*/
	}

	/**
	 * @brief set the possible axis of ration of the body
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	void set_body_rotations(Simulation* simulation, int body_index, bool X, bool Y, bool Z) {
		btVector3 linFact (X ? 1.0f : 0.0f, Y ? 1.0f : 0.0f, Z ? 1.0f : 0.0f);
        simulation->rigidbodies[body_index]->setAngularFactor(linFact);
	}

	/**
	 * @brief Body will move but not rotate (main player)
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	void set_body_no_rotation(Simulation* simulation, int body_index) {
		btVector3 linFact (0.0, 0.0, 0.0);
        simulation->rigidbodies[body_index]->setAngularFactor(linFact);
	}

	glm::vec3 get_body_velocity(Simulation* simulation, int body_index) {
		return bulletToGlm(simulation->rigidbodies[body_index]->getLinearVelocity());
	}

	/**
	 * @brief Set the velocity of the body
	 * 
	 * @param simulation 
	 * @param body_index 
	 * @param velocity 
	 */
	void set_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity) {
		simulation->rigidbodies[body_index]->activate(true);
        simulation->rigidbodies[body_index]->setLinearVelocity(glmToBullet(velocity));
	}

	void set_body_frixion(Simulation* simulation, int body, float frixion) {
		simulation->rigidbodies[body]->setFriction((btScalar)frixion);
	}

	float get_body_frixion(Simulation* simulation, int body) {
		return (float) simulation->rigidbodies[body]->getFriction();
	}

	void set_body_damping(Simulation* simulation, int body, float linear, float angular) {
		simulation->rigidbodies[body]->setDamping((btScalar)linear, (btScalar)angular);
	}

	float get_body_lin_damping(Simulation* simulation, int body) {
		return (float)simulation->rigidbodies[body]->getLinearDamping();
	}

	/**
	 * @brief Get the current velocity of the body and keep only the y component (vertical).
	 * Then it adds the specified velocity to it.
	 * 
	 * @param simulation 
	 * @param body_index 
	 * @param velocity 
	 */
	void add_body_velocity(Simulation* simulation, int body_index, glm::vec3 velocity) {
		simulation->rigidbodies[body_index]->activate(true);
		btVector3 current_vel = simulation->rigidbodies[body_index]->getLinearVelocity();
		current_vel.setX(0.0f);
		current_vel.setZ(0.0f);
		simulation->rigidbodies[body_index]->setLinearVelocity(current_vel + glmToBullet(velocity));
	}

	void print_resume(const Simulation* simulation) {
		std::cout << "Lustrine::Bullet " << "\n" 
			<< "\tnum registered bodies: " << simulation->num_bodies << "\n"
			<< "\tnum registered shapes: " << simulation->collisionShapes.size() << std::endl;
	}

	/**
	 * @brief Apply an impulse on the center of mass of the body
	 * 
	 * @param simulation 
	 * @param body_index 
	 * @param force 
	 */
	void apply_impulse(Simulation* simulation, int body_index, glm::vec3 impulse, glm::vec3 position) {
		simulation->rigidbodies[body_index]->applyImpulse(glmToBullet(impulse), glmToBullet(position));
	}

	/**
	 * @brief helper to print the collisions
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	void print_collisions(Simulation* simulation) {
		for (int i = 0; i < simulation->bodies_collisions.size(); i++) {
			std::cout << i << " -> " << std::endl;
			for (int j = 0; j < simulation->bodies_collisions[i].size(); j++) {
				std::cout << simulation->bodies_collisions[i][j] << ", ";
			}
			std::cout << std::endl;
		}
	}
	/**
	 * @brief store in indices the indices of body colliding with body. Warning, the
	 * pointer is assume to have at least num_bodies entries. Mind the "s".
	 * 
	 * @param simulation 
	 * @param body 
	 * @param indices 
	 * @param size 
	 */
	void check_collisions(Simulation* simulation, int body, int* indices, int* size) {
		memcpy(indices, simulation->bodies_collisions[body].data(), simulation->bodies_collisions[body].size() * 4);
		*size = simulation->bodies_collisions[body].size();
	}

	/**
	 * @brief Get the num bodies object
	 * 
	 * @param simulation 
	 * @return int 
	 */
	int get_num_bodies(Simulation* simulation) {
		return simulation->num_bodies;
	}

	/**
	 * @brief WARNING might not be working need to test
	 * 
	 * @param simulation 
	 * @param body 
	 * @return true 
	 * @return false 
	 */
	bool do_collide(Simulation* simulation, int body) {
		if (simulation->bodies_collisions[body].size() > 0) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @brief returns true if the body considered do collides (are touching in the current frame), false otw.
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
	bool check_collision(Simulation* simulation, int body1, int body2) {
		for (int i = 0; i < simulation->bodies_collisions[body1].size(); i++) {
			if (simulation->bodies_collisions[body1][i] == body2) {
				return true;
			}
		}
		return false;
	}

	/**
	 * @brief Clear and fill the collisions for all the bodies that collided (this function is called in the simulate bullet)
	 * 
	 * @param simulation 
	 * @param body_index 
	 */
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

	/**
	 * @brief returns the current body positions as specified by its motion state.
	 * 
	 * @param simulation 
	 * @param body 
	 * @return glm::vec3 
	 */
	glm::vec3 get_body_position(Simulation* simulation, int body) {
		btTransform t;
		simulation->rigidbodies[body]->getMotionState()->getWorldTransform(t);
		//printf("{%f, %f, %f}\n", t.getOrigin().getX(), t.getOrigin().getY(), t.getOrigin().getZ());
		return bulletToGlm(t.getOrigin());
	}

	/**
	 * @brief Will override the motion state of the object and set its position to be the new
	 * one specified
	 * 
	 * @param simulation 
	 * @param body 
	 * @param new_position 
	 */
	void set_body_position(Simulation* simulation, int body, glm::vec3 new_position) {
		btTransform t;
		simulation->rigidbodies[body]->getMotionState()->getWorldTransform(t);
		t.setOrigin(glmToBullet(new_position));
		simulation->rigidbodies[body]->setWorldTransform(t);
		//simulation->rigidbodies[body]->getMotionState()->setWorldTransform(t);
	}

	/**
	 * @brief Set the particles box colliders positions object (those that are allocated)
	 * 
	 * @param simulation 
	 * @param particles 
	 * @param start 
	 * @param end 
	 */
	void set_particles_box_colliders_positions(Simulation* simulation, glm::vec3* particles, int start_ptr, int end_ptr) {
		
		int num_close = 0;
		for (int i = start_ptr; i < end_ptr; i++) {
			glm::vec3 tmp = particles[i] - simulation->player_position;
			if (glm::distance(particles[i], simulation->player_position) < simulation->player_box_radius) {// (glm::dot(tmp, tmp) < simulation->player_box_radius * simulation->player_box_radius) {
				simulation->rigidbodies[num_close]->setActivationState(ACTIVE_TAG);
				btTransform t;
				simulation->rigidbodies[num_close]->getMotionState()->getWorldTransform(t);
				t.setOrigin(glmToBullet(particles[i]));
				simulation->rigidbodies[num_close]->getMotionState()->setWorldTransform(t);
				simulation->rigidbodies[num_close]->setWorldTransform(t);
				simulation->rigidbodies[num_close]->setInterpolationWorldTransform(t);

				//simulation->rigidbodies[num_close]->setCollisionFlags(simulation->rigidbodies[num_close]->getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
				simulation->rigidbodies[num_close]->setCollisionFlags(simulation->rigidbodies[num_close]->getCollisionFlags() & ~btCollisionObject::CF_NO_CONTACT_RESPONSE);
				simulation->rigidbodies[num_close]->clearForces();
				btVector3 v(0.0, 0.0, 0.0);
				simulation->rigidbodies[num_close]->setInterpolationLinearVelocity(v);
				simulation->rigidbodies[num_close]->setLinearVelocity(v);
				num_close++;
				if (num_close >= simulation->num_particles_allocated) {
					break;
				}
			}
		}
		//std::cout << simulation->player_position.x << std::endl;
		//std::cout << num_close << std::endl;
		for (int i = num_close; i < simulation->num_particles_allocated; i++) {
			simulation->rigidbodies[i]->setCollisionFlags(simulation->rigidbodies[i]->getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
		}

	}


	/**
	 * @brief disable every allocated particles bounding boxes. Meaning that the
	 * bodies will now stop colliding with the sand.
	 * @note WARNING: Because of the interpolation, the boxes could have infinite speed
	 * during a single frame which could cause colliding object (player) to move very fast
	 * as in a huge kinematic collision.
	 * @param simulation 
	 */
	void disable_particles_bounding_boxes(Simulation* simulation) {
		std::cout << "moving out particles bounding boxes" << std::endl;
		for (int i = simulation->ptr_bounding_box_start; i < simulation->ptr_bounding_box_end; i++) {
			//btTransform& t = simulation->sand_particles_colliders[i]->getWorldTransform();
			//t.setOrigin(btVector3(-100.0f, -100.0f, -100.0f)); //TODO HACKY use
        	//simulation->sand_particles_colliders[i]->getMotionState()->setWorldTransform(t);
			//simulation->sand_particles_colliders[i]->setActivationState(DISABLE_SIMULATION);
			simulation->rigidbodies[i]->setCollisionFlags(simulation->rigidbodies[i]->getCollisionFlags() | btCollisionObject::CF_NO_CONTACT_RESPONSE);
		}
	}

	void enable_particles_bounding_boxes(Simulation* simulation) {
		for (int i = simulation->ptr_bounding_box_start; i < simulation->ptr_bounding_box_end; i++) {
			//btTransform& t = simulation->sand_particles_colliders[i]->getWorldTransform();
			//t.setOrigin(btVector3(-100.0f, -100.0f, -100.0f)); //TODO HACKY use
        	//simulation->sand_particles_colliders[i]->getMotionState()->setWorldTransform(t);
			//simulation->sand_particles_colliders[i]->setActivationState(DISABLE_SIMULATION);
			simulation->rigidbodies[i]->setCollisionFlags(simulation->rigidbodies[i]->getCollisionFlags() & ~btCollisionObject::CF_NO_CONTACT_RESPONSE);
		}
	}

	/**
	 * @brief Will bind the positions used interally in the bullet for bounding boxes
	 * with the actual positions from the particle simulation.
	 * 
	 * @param simulation 
	 * @param foreign_positions 
	 */
	void bind_foreign_sand_positions(Simulation* simulation, glm::vec3* foreign_positions) {
		simulation->foreign_sand_positions = foreign_positions;
	}

	/**
	 * @brief NOT USED CURRENTLY
	 * 
	 * @param simulation 
	 * @param body 
	 * @param particles 
	 * @param start 
	 * @param end 
	 * @param particleRadius 
	 */
	/*
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


	}*/
}
}