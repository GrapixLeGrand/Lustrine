#include "LustrineWrapper.hpp"

#include <iostream>
#include <vector>
#include <cstring>

namespace Lustrine {
namespace Wrapper {

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

	//Simulation* simulation;
	void init_simulation_wrapper(const SimulationParameters* parameters, const Grid* wrapped_grids, const Position* wrapped_positions, int num_grids) {
		
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

		std::cout << "end init" << std::endl;

	}

	void init_grid_box_wrapper(const SimulationParameters* parameters, Grid* wrapped, int X, int Y, int Z, int type, Color color) {
		
		std::cout << "init called" << std::endl;
		Lustrine::Grid original_grid;
		Lustrine::init_grid_box(parameters, &original_grid, X, Y, Z, (Lustrine::MaterialType)type, wrapper_to_glm(color));
		grid_to_grid_wrapper(&original_grid, wrapped);

	}

	void simulate(float dt) {
		std::cout << "dt is " << dt << std::endl;
		Lustrine::simulate(simulation, dt);
	}

	void allocate_grid(Grid* grid, int X, int Y, int Z, bool has_per_cell_color) {
		
		std::cout << "allocating memory" << std::endl;
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
			std::memset(grid->colors, 0, 4 * 4 * grid->num_grid_cells);
		}

		grid->color = { 0, 0, 0, 0 };
		grid->type = 0;
		
		grid->num_occupied_grid_cells = 0;
		
	}

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
		
	}

	void cleanup_simulation() {
		std::cout << "clean up simulation" << std::endl;
		delete simulation;
	}

	void say_hello() {
		std::cout << "hello from cpp!" << std::endl;
	}
}
}