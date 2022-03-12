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


	void allocate_simulation_data(SimulationData** data) {
		*data = new SimulationData();
		(*data)->start_dynamic = 0xBEEF;
	}

	void allocate_simulation_parameters(SimulationParameters** parameters) {
		*parameters = new SimulationParameters();
	}

	void free_simulation_data(SimulationData* data) {
		delete data;
	}

	void free_simulation_parameters(SimulationParameters* parameters) {
		delete parameters;
	}

	glm::vec4 wrapper_to_glm(const Color& color) {
		glm::vec4 result (0.0);
		result.r = color.r;
		result.g = color.g;
		result.b = color.b;
		result.a = color.a;
		return result;
	}

	glm::vec3 wrapper_to_glm(const Position& position) {
		glm::vec3 result(0.0);
		result.x = position.x;
		result.y = position.y;
		result.z = position.z;
		return result;
	}

	Position glm_to_wrapper(const glm::vec3& position) {
		Position result{ 0.0f, 0.0f, 0.0f };
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

	void allocate_grids(Grid** grids, int num_grids) {
		__debug_msg("allocated an array of grids");
		*grids = new Grid[num_grids];
	}

	//Simulation* simulation;
	void init_simulation(const SimulationParameters* parameters, SimulationData* data, const Grid* wrapped_grids, const Position* wrapped_positions, int num_grids) {
		
		std::cout << "wrapper init called!" << std::endl;
		
		if (num_grids <= 0) {
			return;
		}

		simulation = new Simulation();

		std::vector<Lustrine::Grid> original_grids (num_grids);
		std::vector<glm::vec3> original_positions(num_grids);
		
		//init the converted grids
		for (int i = 0; i < num_grids; i++) {
			grid_wrapper_to_grid(&wrapped_grids[i], &original_grids[i]);
		}
		
		//convert the positions
		for (int i = 0; i < num_grids; i++) {
			original_positions[i] = wrapper_to_glm(wrapped_positions[i]);
		}

		init_simulation(parameters, simulation, original_grids, original_positions);

		data->num_particles = simulation->num_particles;
		data->start_dynamic = simulation->ptr_fluid_start;
		data->end_dynamic = simulation->ptr_fluid_end;
		
		data->start_static = simulation->ptr_static_start;
		data->end_static = simulation->ptr_static_end;

		std::cout << "end init" << std::endl;

	}

	void init_grid_box(const SimulationParameters* parameters, Grid* wrapped, int X, int Y, int Z, int type, Color color) {
		
		std::cout << "init called" << std::endl;
		Lustrine::Grid original_grid;
		Lustrine::init_grid_box(parameters, &original_grid, X, Y, Z, (Lustrine::MaterialType)type, wrapper_to_glm(color));
		grid_to_grid_wrapper(&original_grid, wrapped);

	}

	void simulate(float dt) {
		Lustrine::simulate(simulation, dt);
	}

	/*
	void allocate_grid(Grid* grid, int X, int Y, int Z, bool has_per_cell_color) {
		__debug_msg("allocating a grid");
		grid->X = X;
		grid->Y = Y;
		grid->Z = Z;

		grid->num_grid_cells = X * Y * Z;

		if (grid->num_grid_cells <= 0) {
			return;
		}

		grid->has_one_color_per_cell = has_per_cell_color;
		grid->cells = new bool[grid->num_grid_cells];
		std::memset(grid->cells, 0, grid->num_grid_cells);
		
		if (grid->has_one_color_per_cell == true) {
			grid->colors = new Color[grid->num_grid_cells];
			std::memset(grid->colors, 0, sizeof(Color) * grid->num_grid_cells);
		}

		grid->color = { 0.0f, 0.0f, 0.0f, 0.0f };
		grid->type = 0;
		
		grid->num_occupied_grid_cells = 0;
		
	}*/
	/*
	void free_grid(Grid* grid) {

		std::cout << "freeing memory" << std::endl;
		if (grid->num_grid_cells <= 0) {
			return;
		}
		
		if (grid->cells != nullptr) {
			delete grid->cells;
			grid->cells = nullptr;
		}
		if (grid->has_one_color_per_cell == true && grid->colors != nullptr) {
			delete grid->colors;
			grid->colors = nullptr;
		}
		
	}*/

	void free_grids(Grid* grids) {
		std::cout << "freeing memory" << std::endl;
		delete grids;
	}

	void cleanup_simulation() {
		std::cout << "clean up simulation" << std::endl;
		delete simulation;
	}

	void simulation_bind_positions(float** position_ptr, int num_positions) {
		*position_ptr = (float*) simulation->positions.data();
	}

	void say_hello() {
		std::cout << "hello from cpp!" << std::endl;
	}

	void simulation_bind_positions_copy(float* position_ptr) {
        std::memcpy(position_ptr, simulation->positions.data(), sizeof(float) * 3 * simulation->num_particles);
	}

    void read_vox_scene(BindingString *data, const uint8_t *buffer, uint32_t size) {
        std::string s{read_vox_scene_json(buffer, 0)};
        data->length = s.size();
        data->data = new char[s.size() + 1];
        memcpy(data->data, s.c_str(), s.size());
        data->data[s.size()] = 0;
	}

	void free_string(BindingString* data) {
	    delete data->data;
	    data->length = 0;
	    data->data = nullptr;
	}
}
}