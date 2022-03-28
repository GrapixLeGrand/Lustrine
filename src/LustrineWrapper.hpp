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

		int num_sand_particles;
		int num_solid_particles;
		
		int start_sand_index;
		int end_sand_index;

		int start_solid_index;
		int end_solid_index;

	};

	struct Color {
		float r;
		float g;
		float b;
		float a;
	};

	struct Vec3 {
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
	extern "C" LUSTRINE_WRAPPER_EXPORT void init_simulation(
		const SimulationParameters* parameters,
		SimulationData* data,
		const Grid* sand_grids, 
		const Vec3* sand_grids_positions, 
		int num_sand_grids,
		const Grid* solid_grids, 
		const Vec3* solid_grids_positions,
		int num_solid_grids
	);

	extern "C" LUSTRINE_WRAPPER_EXPORT void simulate(float dt);
	extern "C" LUSTRINE_WRAPPER_EXPORT void simulation_bind_positions_copy(float* position_ptr);
	extern "C" LUSTRINE_WRAPPER_EXPORT void cleanup_simulation();
	extern "C" LUSTRINE_WRAPPER_EXPORT void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, int type, Color color);
	extern "C" LUSTRINE_WRAPPER_EXPORT void read_vox_scene(BindingString *data, const uint8_t *buffer, int64_t size);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_string(BindingString* data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void create_grid(Grid* grid, const wchar_t* path, int type, int pathlen);

	//physics
	/**
	 * @brief add a box and returns its identifier
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT int add_box(Vec3 position, bool is_dynamic, Vec3 half_dimensions);
	
	/**
	 * @brief returns true if two bodies pointed by indices are colliding
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT bool check_collision(int body1, int body2);

	/**
	 * @brief returns true if it collides with anything
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT bool do_collide(int body);

	/**
	 * @brief write all collisions indices in the given array, Warning, indices must be at least
	 * of num_bodies.
	 * @param body the given body
	 * @param indices a pointer array of indices
	 * @param the amount of recorded collisions
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT void check_collisions(int body, int* indices, int* size);


	extern "C" LUSTRINE_WRAPPER_EXPORT int get_num_bodies();

	/**
	 * @brief apply an impulse to the designated body
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT void apply_impulse(int body, Vec3 impulse, Vec3 relative_pos);
	
	/**
	 * @brief returns the position of the body
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT Vec3 get_position(int body);

	/**
	 * @brief set the velocity of the body
	 * 
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT void set_velocity(int body, Vec3 velocity);

	/**
	 * @brief remove body's rotation
	 * 
	 */
	extern "C" LUSTRINE_WRAPPER_EXPORT void set_body_no_rotation(int body);

	//extern "C" LUSTRINE_WRAPPER_EXPORT void simulation_bind_positions(float** position_ptr, int num_positions);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void say_hello();

	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_data(SimulationData** data);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_parameters(SimulationParameters** parameters);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grids(Grid** grid, int num_grids);

	/*
	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_data(SimulationData** data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_simulation_parameters(SimulationParameters** params);
	extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grids(Grid**, int num_grids);

	extern "C" LUSTRINE_WRAPPER_EXPORT void free_simulation_data(SimulationData* data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_simulation_parameters(SimulationParameters* parameters);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_grids(Grid* grids);
	*/
	
	
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_grid(Grid* grid, int X, int Y, int Z, bool has_per_cell_color);
	//extern "C" LUSTRINE_WRAPPER_EXPORT void free_grid(Grid* grid);

	

} //Wrapper
} //Lustrine
