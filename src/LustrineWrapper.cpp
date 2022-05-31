#include "LustrineWrapper.hpp"
#include "VoxelLoader.hpp"

#include <iostream>
#include <sstream>
#include "../thirdparty/opengametools/ogt_vox.h"
#include <vector>
#include <cstring>
#include <filesystem>
#include <string>
#include "Simulate.hpp"

namespace Lustrine {
namespace Wrapper {

	void __debug_msg(const std::string& msg) {
		std::cout << "Lustrine::Debug" << " " << msg << std::endl;
	}

	/**
	 * @brief Helper method that convert structs from/to glm
	 * 
	 * @param color 
	 * @return glm::vec4 
	 */
	glm::vec4 wrapper_to_glm(const Color& color) {
		glm::vec4 result (0.0);
		result.r = color.r;
		result.g = color.g;
		result.b = color.b;
		result.a = color.a;
		return result;
	}

	/**
	 * @brief Helper method that convert structs from/to glm
	 * 
	 * @param color 
	 * @return glm::vec4 
	 */
	glm::vec3 wrapper_to_glm(const Vec3& position) {
		glm::vec3 result(0.0);
		result.x = position.x;
		result.y = position.y;
		result.z = position.z;
		return result;
	}

	/**
	 * @brief Helper method that convert structs from/to glm
	 * 
	 * @param color 
	 * @return glm::vec4 
	 */
	Vec3 glm_to_wrapper(const glm::vec3& position) {
		Vec3 result{ 0.0f, 0.0f, 0.0f };
		result.x = position.x;
		result.y = position.y;
		result.z = position.z;
		return result;
	}

	/**
	 * @brief Helper method that convert structs from/to glm
	 * 
	 * @param color 
	 * @return glm::vec4 
	 */
	Color glm_to_wrapper(const glm::vec4& color) {
		Color result{ 0.0f, 0.0f, 0.0f, 0.0f };
		result.r = color.r;
		result.g = color.g;
		result.b = color.b;
		result.a = color.a;
		return result;
	}

	/**
	 * @brief Helper method that convert grids from simulation into wrapper and vice versa
	 * The one from Lustrine and the one from LustrineWrapper.
	 * @param color 
	 * @return glm::vec4 
	 */
	void grid_wrapper_to_grid(const GridWrapper* wrapped, Lustrine::Grid* original) {

		original->cells = std::vector<int>(wrapped->num_grid_cells);
		original->has_one_color_per_cell = wrapped->has_one_color_per_cell;
		original->X = wrapped->X;
		original->Y = wrapped->Y;
		original->Z = wrapped->Z;
		original->type = (Lustrine::MaterialType) wrapped->type;
		original->color = wrapper_to_glm(wrapped->color);
		original->num_occupied_grid_cells = wrapped->num_occupied_grid_cells;
		original->num_grid_cells = wrapped->num_grid_cells;
		original->position = wrapper_to_glm(wrapped->position);

		for (int i = 0; i < original->num_grid_cells; i++) {
			original->cells[i] = wrapped->cells[i];
		}

		if (original->has_one_color_per_cell) {
			original->colors = std::vector<glm::vec4>(original->num_grid_cells);
			for (int i = 0; i < original->num_grid_cells; i++) {
				original->colors[i] = wrapper_to_glm(wrapped->colors[i]);
			}
		}

	}

	/**
	 * @brief Helper method that convert grids from simulation into wrapper and vice versa
	 * The one from Lustrine and the one from LustrineWrapper.
	 * @param color 
	 * @return glm::vec4 
	 */
	void grid_to_grid_wrapper(const Lustrine::Grid* original, GridWrapper* wrapped) {

		wrapped->cells = new int[original->num_grid_cells];
		wrapped->has_one_color_per_cell = original->has_one_color_per_cell;
		wrapped->X = original->X;
		wrapped->Y = original->Y;
		wrapped->Z = original->Z;
		wrapped->type = original->type;
		wrapped->color = glm_to_wrapper(original->color);
		wrapped->num_occupied_grid_cells = original->num_occupied_grid_cells;
		wrapped->num_grid_cells = original->num_grid_cells;
		wrapped->position = glm_to_wrapper(original->position);
		for (int i = 0; i < wrapped->num_grid_cells; i++) {
			wrapped->cells[i] = original->cells[i];
		}

		if (wrapped->has_one_color_per_cell) {
			wrapped->colors = new Color[wrapped->num_grid_cells];
			for (int i = 0; i < original->num_grid_cells; i++) {
				wrapped->colors[i] = glm_to_wrapper(original->colors[i]);
			}
		}

	}


	int get_num_sand_particles() {
		return simulation->num_sand_particles;
	}

	void init_grid_magikavoxel(GridWrapper* grid, const char* path, Vec3 position) {
		std::cout << "LustrineWrapper: init grid with magika at " << path << std::endl;
		Lustrine::Grid tmp;
		Lustrine::init_grid_from_magika_voxel(&tmp, path, wrapper_to_glm(position), Lustrine::MaterialType::SOLID);
		grid_to_grid_wrapper(&tmp, grid);
	}

	/**
	 * @brief Init the wrapped simulation
	 * 
	 * @param parameters 
	 * @param data 
	 * @param sand_grids 
	 * @param sand_grids_positions 
	 * @param num_sand_grids 
	 * @param solid_grids 
	 * @param solid_grids_positions 
	 * @param num_solid_grids 
	 */
	void init_simulation(
		const SimulationParameters* parameters,
		SimulationData* data,
		GridWrapper* sand_grid,
		int num_sand_grids,
		GridWrapper* solid_grids, 
		int num_solid_grids,
		int subdivision
	) {
		
		std::cout << "Lustrine::wrapper init called!" << std::endl;
		
		if (num_sand_grids < 0 || num_solid_grids < 0) {
			std::cout << "Lustrine bad arguments...?" << std::endl;
			return;
		}

		//saved_grids = new GridWrapper*[max_num_grids];
		memset(saved_grids, 0, max_num_grids * sizeof(GridWrapper*));
		simulation = new Simulation();
		std::vector<Lustrine::Grid> original_sand_grids(num_sand_grids);
		
		//init the converted grids
		for (int i = 0; i < num_sand_grids; i++) {
			grid_wrapper_to_grid(&sand_grid[i], &original_sand_grids[i]);
			saved_grids[current_saved_grids_num++] = &sand_grid[i];
		}

		//solid
		std::vector<Lustrine::Grid> original_solid_grids(num_solid_grids);
		
		//init the converted grids
		for (int i = 0; i < num_solid_grids; i++) {
			grid_wrapper_to_grid(&solid_grids[i], &original_solid_grids[i]);
			saved_grids[current_saved_grids_num++] = &solid_grids[i];
		}

		Lustrine::init_simulation(
			parameters,
			simulation,
			original_sand_grids,
			original_solid_grids,
			subdivision
		);

		data->start_sand_index = simulation->ptr_sand_start;
		data->end_sand_index = simulation->ptr_sand_end;

		data->start_solid_index = simulation->ptr_solid_start;
		data->end_solid_index = simulation->ptr_solid_end;

		data->num_sand_particles = simulation->num_sand_particles;
		data->num_solid_particles = simulation->num_solid_particles;

		std::cout << "Lustrine::wrapper: end init" << std::endl;

	}

	void init_simulation_extra_parameters(
		const SimulationParameters* parameters, 
		SimulationData* data, 
		GridWrapper* sand_grids, 
		int num_sand_grids, 
		GridWrapper* solid_grids, 
		int num_solid_grids, 
		int subdivision, 
		float kernel_radius_scale,
		int with_credits
	) {
		std::cout << "Lustrine::wrapper init called!" << std::endl;

		if (num_sand_grids < 0 || num_solid_grids < 0) {
			std::cout << "Lustrine bad arguments...?" << std::endl;
			return;
		}

		//saved_grids = new GridWrapper*[max_num_grids];
		memset(saved_grids, 0, max_num_grids * sizeof(GridWrapper*));
		simulation = new Simulation();
		std::vector<Lustrine::Grid> original_sand_grids(num_sand_grids);

		//init the converted grids
		for (int i = 0; i < num_sand_grids; i++) {
			grid_wrapper_to_grid(&sand_grids[i], &original_sand_grids[i]);
			saved_grids[current_saved_grids_num++] = &sand_grids[i];
		}

		//solid
		std::vector<Lustrine::Grid> original_solid_grids(num_solid_grids);

		//init the converted grids
		for (int i = 0; i < num_solid_grids; i++) {
			grid_wrapper_to_grid(&solid_grids[i], &original_solid_grids[i]);
			saved_grids[current_saved_grids_num++] = &solid_grids[i];
		}

		Lustrine::init_simulation_extra_parameters(
			parameters,
			simulation,
			original_sand_grids,
			original_solid_grids,
			subdivision,
			kernel_radius_scale,
			with_credits
		);

		data->start_sand_index = simulation->ptr_sand_start;
		data->end_sand_index = simulation->ptr_sand_end;

		data->start_solid_index = simulation->ptr_solid_start;
		data->end_solid_index = simulation->ptr_solid_end;

		data->num_sand_particles = simulation->num_sand_particles;
		data->num_solid_particles = simulation->num_solid_particles;

		std::cout << "Lustrine::wrapper: end init" << std::endl;

	}

	void set_player_particles_bounding_spheres_radius_placement(float radius) {
		simulation->bullet_physics_simulation.player_box_radius = radius;
	}

	int query_cell_num_particles(Vec3 min, Vec3 max, bool include_solid)
	{
		return Lustrine::query_cell_num_particles(simulation, wrapper_to_glm(min), wrapper_to_glm(max), include_solid);
	}

	/**
	 * @brief Get the current gravity of the bullet world
	 * 
	 * @return Vec3 
	 */
	Vec3 get_gravity() {
		return glm_to_wrapper(Lustrine::Bullet::get_gravity(&simulation->bullet_physics_simulation));
	}

	/**
	 * @brief Set the gravity of the bullet physics simulation
	 * 
	 * @param new_gravity 
	 */
	void set_gravity(Vec3 new_gravity) {
		Lustrine::Bullet::set_gravity(&simulation->bullet_physics_simulation, wrapper_to_glm(new_gravity));
	}

	/**
	 * @brief add a capsule of height and radius.
	 * The capsule is automatically dynamic (will fall).
	 * @return the body index of the new capsule
	 * 
	 */
	int add_capsule(Vec3 position, float radius, float height) {
		return Lustrine::Bullet::add_capsule(&simulation->bullet_physics_simulation, wrapper_to_glm(position), radius, height);
	}

	/**
	 * @brief Adds a detector block at position and of dimensions halfdims * 2
	 * @return the index of the new detector block
	 */
	int add_detector_block(Vec3 position, Vec3 half_dims) {
		int id = Lustrine::Bullet::add_detector_block(&simulation->bullet_physics_simulation, wrapper_to_glm(position), wrapper_to_glm(half_dims));
		set_body_gravity(id, { 0.0f, 0.0f, 0.0f });
		return id;
	}

	//int add_detector_cylinder(Vec3 position, Vec3 half_dims)
	//{
	//	return Lustrine::Bullet::add_detector_cylinder(&simulation->bullet_physics_simulation, wrapper_to_glm(position), wrapper_to_glm(half_dims));
	//}

	void init_grid_box(const SimulationParameters* parameters, GridWrapper* wrapped, int X, int Y, int Z, Vec3 position, Color color, int type) {
		std::cout << "init grid called" << std::endl;
		Lustrine::Grid original_grid;
		Lustrine::init_grid_box(parameters, &original_grid, X, Y, Z, wrapper_to_glm(position), wrapper_to_glm(color), (Lustrine::MaterialType)type);
		grid_to_grid_wrapper(&original_grid, wrapped);
	}

	void simulate(float dt, bool attract_flag, bool blow_flag) {
		simulation->attract_flag = attract_flag;
		simulation->blow_flag = blow_flag;
		Lustrine::simulate(simulation, dt); //, wrapper_to_glm(character_pos), attract_flag, blow_flag);
	}
	
	void simulate_no_flags(float dt) {
		Lustrine::simulate(simulation, dt);
	}
		
	void cleanup_simulation() {
		std::cout << "clean up simulation" << std::endl;
		Lustrine::clean_simulation(simulation);
		delete simulation;
		simulation = nullptr;

		/*
		std::cout << "cleaning up the " << current_saved_grids_num << " Grids\n";
		for (int i = 0; i < current_saved_grids_num; i++) {
			GridWrapper* current = saved_grids[i];
			if (current == nullptr) {
				std::cout << "null grid encountered on cleanup\n";
			}
			std::cout << "here " << std::endl;
			delete[] current->cells;
			if (current->has_one_color_per_cell) {
				std::cout << "color" << std::endl;
				delete[] current->colors;
			}
			
		}*/

	}

	void simulation_bind_positions_copy(float* position_ptr) {
        std::memcpy(position_ptr, simulation->positions, sizeof(float) * 3 * simulation->num_sand_particles);
	}

    void read_vox_scene(BindingString *data, const uint8_t *buffer, int64_t size) {
        std::string json{read_vox_scene_json(buffer, size)};
        data->length = json.size();
        data->data = new char[json.size() + 1];
        memcpy(data->data, json.c_str(), json.size());
        data->data[json.size()] = 0;
	}

	void free_string(BindingString* data) {
	    delete data->data;
	    data->length = 0;
	    data->data = nullptr;
	}

	/**
	 * @brief write all collisions indices in the given array, Warning, indices must be at least
	 * of num_bodies.
	 * @param body the given body
	 * @param indices a pointer array of indices
	 * @param the amount of recorded collisions
	 */
	void check_collisions(int body, int* indices, int* size) {
		Lustrine::Bullet::check_collisions(&simulation->bullet_physics_simulation, body, indices, size);
	}

	int get_num_bodies() {
		return Lustrine::Bullet::get_num_bodies(&simulation->bullet_physics_simulation);
	}

	/**
	 * @brief add a box and returns its identifier
	 */
	int add_box(Vec3 position, bool is_dynamic, Vec3 half_dimensions) {
		return Lustrine::Bullet::add_box(&simulation->bullet_physics_simulation, wrapper_to_glm(position), is_dynamic, wrapper_to_glm(half_dimensions)); //, -1, INT32_MAX);
	}

	/**
	 * @brief returns true if two bodies pointed by indices are colliding
	 */
	int check_collision(int body1, int body2) {
		return (int) Lustrine::Bullet::check_collision(&simulation->bullet_physics_simulation, body1, body2);
	}

	/**
	 * @brief returns true if the body collides with anything
	 */
	int do_collide(int body) {
		return (int) Lustrine::Bullet::do_collide(&simulation->bullet_physics_simulation, body);
	}

	int do_collide_except_for(int body, int exception_id) {
		return (int)Lustrine::Bullet::do_collide_except_for(&simulation->bullet_physics_simulation, body, exception_id);
	}

	/**
	 * @brief apply an impulse to the designated body
	 */
	void apply_impulse(int body, Vec3 impulse, Vec3 relative_pos) {
		Lustrine::Bullet::apply_impulse(&simulation->bullet_physics_simulation, body, wrapper_to_glm(impulse), wrapper_to_glm(relative_pos));
	}

	/**
	 * @brief returns the position of the body
	 */
	Vec3 get_position(int body) {
		return glm_to_wrapper(Lustrine::Bullet::get_body_position(&simulation->bullet_physics_simulation, body));
	}

	glm::vec3 get_velocity(int body) {
		return Lustrine::Bullet::get_body_velocity(&simulation->bullet_physics_simulation, body);
	}
	/**
	 * @brief set the velocity of the body
	 * 
	 */
	void set_velocity(int body, Vec3 velocity) {
		Lustrine::Bullet::set_body_velocity(&simulation->bullet_physics_simulation, body, wrapper_to_glm(velocity));
	}

	void add_velocity(int body, Vec3 velocity) {
		Lustrine::Bullet::add_body_velocity(&simulation->bullet_physics_simulation, body, wrapper_to_glm(velocity));
	}

	void set_position(int body, Vec3 position) {
		Lustrine::Bullet::set_body_position(&simulation->bullet_physics_simulation, body, wrapper_to_glm(position));
	}

	void set_body_frixion(int body, float frixion) {
		Lustrine::Bullet::set_body_frixion(&simulation->bullet_physics_simulation, body, frixion);
	}

	extern "C" LUSTRINE_WRAPPER_EXPORT void set_body_damping(int body, float linear, float angular) {
		Lustrine::Bullet::set_body_damping(&simulation->bullet_physics_simulation, body, linear, angular);
	}

	extern "C" LUSTRINE_WRAPPER_EXPORT float get_body_damping(int body) {
		return Lustrine::Bullet::get_body_lin_damping(&simulation->bullet_physics_simulation, body);
	}

	/**
	 * @brief set player's bullet id
	 *
	 */
	void set_player_id(int id)
	{
		simulation->bullet_physics_simulation.player_id = id;
	}

	void set_player_box_scale(Vec3 scale)
	{
		simulation->bullet_physics_simulation.player_box_scale = wrapper_to_glm(scale);
	}

	int is_grounded(int id)
	{
		glm::vec3 pos = Bullet::get_body_position(&simulation->bullet_physics_simulation, id);
		btVector3 btFrom(pos.x, pos.y, pos.z);
		btVector3 btTo(pos.x, pos.y - 0.55f, pos.z);
		btCollisionWorld::ClosestRayResultCallback rayCallback(btFrom, btTo);
		simulation->bullet_physics_simulation.dynamicWorld->rayTest(btFrom, btTo, rayCallback);

		return (int)rayCallback.hasHit();
	}

	void set_attract_blow_parameters(float attract_radius, float blow_radius, float attract_coeff, float blow_coeff)
	{
		simulation->attract_radius = attract_radius;
		simulation->blow_radius = blow_radius;
		simulation->attract_coeff = attract_coeff;
		simulation->blow_coeff = blow_coeff;
	}


	void set_body_gravity(int id, Vec3 gravity)
	{
		glm::vec3 gravity_tmp = wrapper_to_glm(gravity);
		simulation->bullet_physics_simulation.rigidbodies[id]->setGravity(btVector3(gravity_tmp.x, gravity_tmp.y, gravity_tmp.z));
	}

	void set_body_no_collision_response(int id)
	{
		simulation->bullet_physics_simulation.rigidbodies[id]->setCollisionFlags(btCollisionObject::CF_NO_CONTACT_RESPONSE);
	}

	int collide_with_player(int id)
	{
		return check_collision(id, simulation->bullet_physics_simulation.player_id);
	}
	/**
	 * @brief remove body's rotation
	 * 
	 */
	void set_body_no_rotation(int body) {
		Lustrine::Bullet::set_body_no_rotation(&simulation->bullet_physics_simulation, body);
	}

	void create_grid(GridWrapper* grid, const wchar_t* path, int type, int pathlen) {
		
		// Rewrite of init_grid_from_magika_voxel
		char* newpath = (char*)malloc(pathlen + 1);
		int i;
		for (i = 0; i < pathlen; ++i) {
			newpath[i] = path[i];
		}
		newpath[i] = 0;
		//std::string s(newpath);
		//printf("\n ***********   paththt is %s \n ",newpath);
		FILE* fp;
		try {
			fp = fopen(newpath, "rb");
			free(newpath);
		}
		catch (int e) {
			std::cout << "An exception occurred" << std::endl;
		}
		//std::cout << sizeof(path[0]) << std::endl;
		if (!fp) {
			std::cout << "File could not be allocated" << std::endl;
			return;
		}
		
		// get the buffer size which matches the size of the file
		fseek(fp, 0, SEEK_END);
		uint32_t buffer_size = ftell(fp);
		fseek(fp, 0, SEEK_SET);

		// load the file into a memory buffer
		uint8_t* buffer = new uint8_t[buffer_size];
		fread(buffer, buffer_size, 1, fp);
		fclose(fp);

		// construct the scene from the buffer
		const ogt_vox_scene* scene = ogt_vox_read_scene_with_flags(buffer, buffer_size, 0);
		
		// the buffer can be safely deleted once the scene is instantiated.
		delete[] buffer;
		
		if (scene->num_models < 1) {
			std::cout << "Emtpy voxel model" << std::endl;
		}
		
		std::cout << scene->num_layers << std::endl;
		const ogt_vox_model* model = scene->models[0];

		grid->X = model->size_x;
		grid->Y = model->size_y;
		grid->Z = model->size_z;

		grid->type = type;
		grid->num_grid_cells = model->size_x * model->size_y * model->size_z;
		grid->cells = (int*)calloc(grid->num_grid_cells,sizeof(bool)); //std::vector<bool>(grid->num_grid_cells, false);
		grid->colors = (Color*)calloc(grid->num_grid_cells, sizeof(Color));  //std::vector<glm::vec4>(grid->num_grid_cells, { 0, 0, 0, 0 });
		grid->has_one_color_per_cell = true;
		
		int counter = 0;
		for (int x = 0; x < grid->X; x++) {
			for (int y = 0; y < grid->Y; y++) {
				for (int z = 0; z < grid->Z; z++) {
					int voxel_index = x + (y * model->size_x) + (z * model->size_x * model->size_y);
					int grid_index = (x * grid->Y * grid->Z) + (y * grid->Z) + z;
					uint8_t color_index = model->voxel_data[voxel_index];

					if (color_index == 0) { //voxel is non existent
						continue;
					}

					ogt_vox_rgba voxel_color = scene->palette.color[color_index];
					grid->cells[grid_index] = true;
					Color& grid_color = grid->colors[grid_index];

					grid_color.r = voxel_color.r/255.0;
					grid_color.g = voxel_color.g/255.0;
					grid_color.b = voxel_color.b/255.0;
					grid_color.a = voxel_color.a/255.0;
					//grid_color /= 255.0;

					counter++;
				}
			}
		}

		grid->num_occupied_grid_cells = counter;
		ogt_vox_destroy_scene(scene);

		return;
	}

	/**
	 * @brief allocate particles collider. Meaning that we allocate a box (for now at least)
	 * that will encapsulate the particles. This functions only allocates the boxes but never
	 * sets the actual positions of the bodies.
	 * @note This function can be recalled with another argument, possibly resizing the existing
	 * allocation.
	 * @param num_particles 
	 */
	/*
	void allocate_particles_colliders() {
		Lustrine::Bullet::allocate_particles_colliders(&simulation->bullet_physics_simulation);
	}*/

	/**
	 * @brief Set the bounding box (in bullet) of the particles to the particles
	 * position (in the simulation)
	 */
	/*
	void set_particles_box_colliders_positions() {
		Lustrine::Bullet::set_particles_box_colliders_positions(&simulation->bullet_physics_simulation, simulation->positions);
	}*/

	/**
	 * @brief enables the bounding box of the particles (can collide with the player)
	 * 
	 */
	void enable_particles_bounding_boxes() {
		simulation->bullet_physics_simulation.particles_bounding_box_requested_state = true;
	}

	/**
	 * @brief disable the bounding boxes of the particles (cannot collide with the player)
	 * 
	 */
	void disable_particles_bounding_boxes() {
		simulation->bullet_physics_simulation.particles_bounding_box_requested_state = false;
	}

	void set_simulate_function(int index) {
		switch (index) {
			case 0:
			simulation->simulate_fun = simulate_sand;
			break;
			case 1:
			simulation->simulate_fun = simulate_sand_v3;
			break;
			default:
			std::cout << "Unreckognized input for simulate func " << index << "\n";
		}
	}

	int add_particle_source(GridWrapper* pattern, Vec3 direction, float freq, int capacity) {
		Lustrine::Grid tmp;
		grid_wrapper_to_grid(pattern, &tmp);
		return Lustrine::add_particle_source(simulation, &tmp, wrapper_to_glm(direction), freq, capacity);
	}

	int add_particle_sink(Vec3 min_pos, Vec3 max_pos, float frequency) {
		return Lustrine::add_particle_sink(simulation, wrapper_to_glm(min_pos), wrapper_to_glm(max_pos), frequency);
	}

	void set_source_state(int index, int state) {
		Lustrine::set_source_state(simulation, index, state != 0 ? true : false);
	}

	void set_sink_state(int index, int state) {
		Lustrine::set_sink_state(simulation, index, state != 0 ? true : false);
	}

	int get_source_spawned(int index) {
		return Lustrine::get_source_spawned(simulation, index);
	}

	int get_sink_despawned(int index) {
		return Lustrine::get_sink_despawned(simulation, index);
	}

	int get_grid_cell_size() {
		return simulation->cell_size;
	}


	static char* dummy_ptr = nullptr;
	std::align_val_t simd_vector_align{ 64 };

	extern "C" LUSTRINE_WRAPPER_EXPORT void test_allocate_1gb() {
		dummy_ptr = new (simd_vector_align) char[1000000000];
	}

	extern "C" LUSTRINE_WRAPPER_EXPORT void test_deallocate_1gb() {
		::operator delete[](dummy_ptr, simd_vector_align);
	}
}	
}