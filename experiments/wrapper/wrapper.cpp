#include <iostream>

#include "LustrineWrapper.hpp"

int main(void) {
	std::cout << "Wrapper experiment" << std::endl;
	Lustrine::SimulationParameters parameters;
	parameters.X = 30;
	parameters.Y = 20;
	parameters.Z = 30;

	Lustrine::Wrapper::SimulationData data;
	Lustrine::Wrapper::Grid grids[2];
	Lustrine::Wrapper::Position positions[2];

	positions[0] = {0, 0, 0};
	positions[1] = { 10, 0, 10 };

	Lustrine::Wrapper::init_grid_box(&parameters, &grids[0], 10, 10, 10, 0, { 1.0f, 0.0f, 0.0f, 1.0f });
	Lustrine::Wrapper::init_grid_box(&parameters, &grids[1], 10, 10, 10, 0, {1.0f, 0.0f, 0.0f, 1.0f});

	Lustrine::Wrapper::init_simulation(&parameters, &data, grids, positions, 2);

	std::cout << "simulate for 1000 steps !" << std::endl;

	float* particles_positions = nullptr;
	
	int p = 0;
	Lustrine::Wrapper::simulation_bind_positions(&particles_positions, p);
	const int size = 3 * data.num_particles;
	float* test = new float[size];
	Lustrine::Wrapper::simulation_bind_positions_copy(test);

	for (int i = 0; i < 2000; i++) {
		std::cout << particles_positions[3 * i + 0] << particles_positions[3 * i + 1] << particles_positions[3 * i + 2]  << std::endl;
	}

	float dt = 0.01;
	for (int i = 0; i < 1000; i++) {
		Lustrine::Wrapper::simulate(dt);
		std::cout << particles_positions[33] << " " << particles_positions[34] << " " << particles_positions[35] << std::endl;
	}

	std::cout << "dt = " << dt << std::endl;
	Lustrine::Wrapper::cleanup_simulation();

}