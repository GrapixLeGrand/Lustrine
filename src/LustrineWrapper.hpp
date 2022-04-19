#pragma once

/**
 * @file LustrineWrapper.hpp
 * @author Quentin Guignard (qguignard@student.ethz.ch), Gilles Waeber
 * @brief 
 * @version 0.1
 * @date 2022-04-06
 * 
 * @copyright Copyright (c) 2022
 * 
 */

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

	/**
	 * @brief Wrapper structure for Lustrine::Grid. WARNING: TWO grids are defined in the project.
	 * This one in LustrineWrapper (LustrineWrapper::Grid) that wraps the other Lustrine::Grid in
	 * Simulation.hpp. This Grid has NO managed memory and MUST never have some to stay compatible
	 * with the binding.
	 */
	struct Grid {
		bool* cells;
		Color* colors;
		Color color;
		Vec3 position;
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
		int num_sand_grids,
		const Grid* solid_grids,
		int num_solid_grids
	);

	/**
	 * @brief Utils function of the bindings
	 * 
	 */

	extern "C" LUSTRINE_WRAPPER_EXPORT void simulate(float dt);
	extern "C" LUSTRINE_WRAPPER_EXPORT void simulation_bind_positions_copy(float* position_ptr);
	extern "C" LUSTRINE_WRAPPER_EXPORT void cleanup_simulation();
	extern "C" LUSTRINE_WRAPPER_EXPORT void init_grid_box(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, Vec3 position, Color color, int type);
	extern "C" LUSTRINE_WRAPPER_EXPORT void read_vox_scene(BindingString *data, const uint8_t *buffer, int64_t size);
	extern "C" LUSTRINE_WRAPPER_EXPORT void free_string(BindingString* data);
	extern "C" LUSTRINE_WRAPPER_EXPORT void create_grid(Grid* grid, const wchar_t* path, int type, int pathlen);
	extern "C" LUSTRINE_WRAPPER_EXPORT void init_grid_magikavoxel(Grid* grid, const char* path, Vec3 position);

	/**
	 * @brief bullet physics functions. doc in cpp. 
	 * 
	 */
	
	extern "C" LUSTRINE_WRAPPER_EXPORT Vec3 get_gravity();
	extern "C" LUSTRINE_WRAPPER_EXPORT void set_gravity(Vec3 new_gravity);
	extern "C" LUSTRINE_WRAPPER_EXPORT int add_box(Vec3 position, bool is_dynamic, Vec3 half_dimensions);
	extern "C" LUSTRINE_WRAPPER_EXPORT int add_capsule(Vec3 position, float radius, float height);
	extern "C" LUSTRINE_WRAPPER_EXPORT int add_detector_block(Vec3 position, Vec3 half_dims);
	extern "C" LUSTRINE_WRAPPER_EXPORT bool check_collision(int body1, int body2);
	extern "C" LUSTRINE_WRAPPER_EXPORT bool do_collide(int body);
	extern "C" LUSTRINE_WRAPPER_EXPORT void check_collisions(int body, int* indices, int* size);
	extern "C" LUSTRINE_WRAPPER_EXPORT int get_num_bodies();
	extern "C" LUSTRINE_WRAPPER_EXPORT void apply_impulse(int body, Vec3 impulse, Vec3 relative_pos);
	extern "C" LUSTRINE_WRAPPER_EXPORT Vec3 get_position(int body);
	extern "C" LUSTRINE_WRAPPER_EXPORT void set_velocity(int body, Vec3 velocity);
	extern "C" LUSTRINE_WRAPPER_EXPORT void set_body_no_rotation(int body);

	/**
	 * @brief For particles interactions
	 */
	//extern "C" LUSTRINE_WRAPPER_EXPORT void allocate_particles_colliders(int num_particles);
	extern "C" LUSTRINE_WRAPPER_EXPORT void enable_particles_bounding_boxes();
	extern "C" LUSTRINE_WRAPPER_EXPORT void disable_particles_bounding_boxes();

} //Wrapper
} //Lustrine
