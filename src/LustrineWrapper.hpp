#pragma once

#include "LustrineWrapper_Export.h"
#include "Lustrine.hpp"

#ifdef LUSTRINEWRAPPER_API
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define LUSTRINEWRAPPER_API __declspec(dllimport)
#endif

namespace Lustrine {
namespace Wrapper {

	/**
	* This struct is meant to communicate to the client the state of the simulation
	* it should be set by the simulation on init not manually! Don't write in it.
	*/
	struct SimulationData {

		int num_particles;
		
		int start_dynamic;
		int end_dynamic;

		int start_static;
		int end_static;

	};

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

	struct BindingString {
	    int64_t length;
	    char* data;
	};

	static Simulation* simulation = nullptr;
	extern "C" LUSTRINE_WRAPPER_EXPORT void init_simulation(const SimulationParameters* parameters, SimulationData* data, const Grid* grids, const Position* positions, int num_grids);
	extern "C" LUSTRINE_WRAPPER_EXPORT void simulate(float dt);
	extern "C" LUSTRINE_WRAPPER_EXPORT void simulation_bind_positions(float** position_ptr, int num_positions);

	extern "C" LUSTRINE_WRAPPER_EXPORT void simulation_bind_positions_copy(float* position_ptr);

	extern "C" LUSTRINE_WRAPPER_EXPORT void cleanup_simulation();
	extern "C" LUSTRINE_WRAPPER_EXPORT void say_hello();

	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_data(SimulationData** data);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_parameters(SimulationParameters** parameters);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grids(Grid** grid, int num_grids);

	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_data(SimulationData** data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_parameters(SimulationParameters** params);
	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grids(Grid**, int num_grids);

	extern "C" LUSTRINE_WRAPPER_EXPORT void free_simulation_data(SimulationData* data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_simulation_parameters(SimulationParameters* parameters);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_grids(Grid* grids);

	extern "C" LUSTRINE_WRAPPER_EXPORT void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, int type, Color color);
	
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grid(Grid* grid, int X, int Y, int Z, bool has_per_cell_color);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void free_grid(Grid* grid);

	extern "C" LUSTRINE_WRAPPER_EXPORT void read_vox_scene(BindingString *data, const uint8_t *buffer, int64_t size);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_string(BindingString* data);

} //Wrapper
} //Lustrine
