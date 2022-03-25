#include "LustrineWrapper.hpp"
#include "VoxelLoader.hpp"

#include <iostream>
#include <vector>
#include <cstring>
#include <string>

namespace Lustrine {
namespace Wrapper {

	void __debug_msg(const std::string& msg) {
		std::cout << "Lustrine::Debug" << " " << msg << std::endl;
	}

	/*
	void allocate_simulation_data(SimulationData** data) {
		*data = new SimulationData();
		(*data)->start_sand_index = 0xBEEF;
	}

	void allocate_simulation_parameters(SimulationParameters** parameters) {
		*parameters = new SimulationParameters();
	}

	void free_simulation_data(SimulationData* data) {
		delete data;
	}
	
	void free_simulation_parameters(SimulationParameters* parameters) {
		delete parameters;
	}*/

	glm::vec4 wrapper_to_glm(const Color& color) {
		glm::vec4 result (0.0);
		result.r = color.r;
		result.g = color.g;
		result.b = color.b;
		result.a = color.a;
		return result;
	}

	glm::vec3 wrapper_to_glm(const Vec3& position) {
		glm::vec3 result(0.0);
		result.x = position.x;
		result.y = position.y;
		result.z = position.z;
		return result;
	}

	Vec3 glm_to_wrapper(const glm::vec3& position) {
		Vec3 result{ 0.0f, 0.0f, 0.0f };
		result.x = position.x;
		result.y = position.y;
		result.z = position.z;
		return result;
	}

	Color glm_to_wrapper(const glm::vec4& color) {
		Color result{ 0.0f, 0.0f, 0.0f, 0.0f };
		result.r = color.r;
		result.g = color.g;
		result.b = color.b;
		result.a = color.a;
		return result;
	}

	void grid_wrapper_to_grid(const Grid* wrapped, Lustrine::Grid* original) {

		original->cells = std::vector<bool>(wrapped->num_grid_cells);
		original->has_one_color_per_cell = wrapped->has_one_color_per_cell;
		original->X = wrapped->X;
		original->Y = wrapped->Y;
		original->Z = wrapped->Z;
		original->type = (Lustrine::MaterialType) wrapped->type;
		original->color = wrapper_to_glm(wrapped->color);
		original->num_occupied_grid_cells = wrapped->num_occupied_grid_cells;
		original->num_grid_cells = wrapped->num_grid_cells;
		
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

	void grid_to_grid_wrapper(const Lustrine::Grid* original, Grid* wrapped) {

		wrapped->cells = new bool[original->num_grid_cells];
		wrapped->has_one_color_per_cell = original->has_one_color_per_cell;
		wrapped->X = original->X;
		wrapped->Y = original->Y;
		wrapped->Z = original->Z;
		wrapped->type = original->type;
		wrapped->color = glm_to_wrapper(original->color);
		wrapped->num_occupied_grid_cells = original->num_occupied_grid_cells;
		wrapped->num_grid_cells = original->num_grid_cells;

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

	/*
	void allocate_grids(Grid** grids, int num_grids) {
		__debug_msg("allocated an array of grids");
		*grids = new Grid[num_grids];
	}*/

	void init_simulation(
		const SimulationParameters* parameters,
		SimulationData* data,
		const Grid* sand_grids, 
		const Vec3* sand_grids_positions, 
		int num_sand_grids,
		const Grid* solid_grids, 
		const Vec3* solid_grids_positions,
		int num_solid_grids
	) {
		
		std::cout << "Lustrine::wrapper init called!" << std::endl;
		
		if (num_sand_grids < 0 || num_solid_grids < 0) {
			std::cout << "Lustrine bad arguments..." << std::endl;
			return;
		}

		simulation = new Simulation();

		std::vector<Lustrine::Grid> original_sand_grids(num_sand_grids);
		std::vector<glm::vec3> original_sand_positions(num_sand_grids);
		
		//init the converted grids
		for (int i = 0; i < num_sand_grids; i++) {
			grid_wrapper_to_grid(&sand_grids[i], &original_sand_grids[i]);
		}
		
		//convert the positions
		for (int i = 0; i < num_sand_grids; i++) {
			original_sand_positions[i] = wrapper_to_glm(sand_grids_positions[i]);
		}

		//solid
		std::vector<Lustrine::Grid> original_solid_grids(num_solid_grids);
		std::vector<glm::vec3> original_solid_positions(num_solid_grids);
		
		//init the converted grids
		for (int i = 0; i < num_solid_grids; i++) {
			grid_wrapper_to_grid(&solid_grids[i], &original_solid_grids[i]);
		}
		
		//convert the positions
		for (int i = 0; i < num_solid_grids; i++) {
			original_solid_positions[i] = wrapper_to_glm(solid_grids_positions[i]);
		}

		Lustrine::init_simulation(
			parameters,
			simulation,
			original_sand_grids,
			original_sand_positions,
			original_solid_grids,
			original_solid_positions
		);

		data->start_sand_index = simulation->ptr_sand_start;
		data->end_sand_index = simulation->ptr_sand_end;

		data->start_solid_index = simulation->ptr_solid_start;
		data->end_solid_index = simulation->ptr_solid_end;

		data->num_sand_particles = simulation->num_sand_particles;
		data->num_solid_particles = simulation->num_solid_particles;

		std::cout << "Lustrine::wrapper: end init" << std::endl;

	}

	void init_grid_box(const SimulationParameters* parameters, Grid* wrapped, int X, int Y, int Z, int type, Color color) {
		std::cout << "init grid called" << std::endl;
		Lustrine::Grid original_grid;
		Lustrine::init_grid_box(parameters, &original_grid, X, Y, Z, (Lustrine::MaterialType)type, wrapper_to_glm(color));
		grid_to_grid_wrapper(&original_grid, wrapped);
	}

	void simulate(float dt) {
		Lustrine::simulate(simulation, dt);
	}

	void cleanup_simulation() {
		std::cout << "clean up simulation" << std::endl;
		Lustrine::clean_simulation(simulation);
		delete simulation;
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


	int add_box(Vec3 position, bool is_dynamic, int lenX, int lenY, int lenZ) {
		glm::vec3 half = glm::vec3{ lenX, lenY, lenZ };
		return Lustrine::Bullet::add_box(&simulation->bullet_physics_simulation, wrapper_to_glm(position), is_dynamic, half, -1, INT32_MAX);
	}

	bool check_collision(int body1, int body2) {
		return Lustrine::Bullet::check_collision(&simulation->bullet_physics_simulation, body1, body2);
	}

	void apply_impulse(int body, Vec3 impulse, Vec3 relative_pos) {
		Lustrine::Bullet::apply_impulse(&simulation->bullet_physics_simulation, body, wrapper_to_glm(impulse), wrapper_to_glm(relative_pos));
	}

	Vec3 get_position(int body) {
		glm_to_wrapper(Lustrine::Bullet::get_body_position(&simulation->bullet_physics_simulation, body));
	}

	void set_velocity(int body, Vec3 velocity) {
		Lustrine::Bullet::set_body_velocity(&simulation->bullet_physics_simulation, body, wrapper_to_glm(velocity));
	}

	void set_body_no_rotation(int body) {
		Lustrine::Bullet::set_body_no_rotation(&simulation->bullet_physics_simulation, body);
	}
	
}	
}