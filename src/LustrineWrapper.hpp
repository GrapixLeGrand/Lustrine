#pragma once

#include "Lustrine.hpp"

namespace Lustrine {
namespace Wrapper {
	/*
	Wrapping
	struct Grid {
		std::vector<bool> cells;
		std::vector<glm::vec4> colors;
		glm::vec4 color;
		bool has_one_color_per_cell;
		int X;
		int Y;
		int Z;
		int num_grid_cells;
		int num_occupied_grid_cells;
		MaterialType type;
	};
	*/

	struct Color {
		float r;
		float g;
		float b;
		float a;
	};

	struct Position {
		float x;
		float y;
		float z;
	};

	struct Grid {
		bool* cells;
		Color* colors;
		Color color;
		bool has_one_color_per_cell;
		int X;
		int Y;
		int Z;
		int num_grid_cells;
		int num_occupied_grid_cells;
		int type;
	};

	static Simulation* simulation;
	extern void init_simulation_wrapper(const SimulationParameters* parameters, Grid* grids, Position* positions);
	extern void init_grid_box_wrapper(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, int type, Color color);
	extern void simulate_wrapper(float dt);
} //Wrapper
} //Lustrine