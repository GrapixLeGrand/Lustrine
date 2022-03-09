#pragma once

#include "Lustrine.hpp"

#pragma once

#ifdef LUSTRINEWRAPPER_API
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define LUSTRINEWRAPPER_API __declspec(dllimport)
#endif

namespace Lustrine {
namespace Wrapper {

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

	static Simulation* simulation = nullptr;
	extern "C" __declspec(dllexport) void init_simulation_wrapper(const SimulationParameters* parameters, const Grid* grids, const Position* positions, int num_grids);
	extern "C" __declspec(dllexport) void init_grid_box_wrapper(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, int type, Color color);
	extern "C" __declspec(dllexport) void simulate(float dt);

	extern "C" __declspec(dllexport) void allocate_grid(Grid* grid, int X, int Y, int Z, bool has_per_cell_color);
	extern "C" __declspec(dllexport) void free_grid(Grid* grid);
	extern "C" __declspec(dllexport) void cleanup_simulation();

	extern "C" __declspec(dllexport) void say_hello();

} //Wrapper
} //Lustrine