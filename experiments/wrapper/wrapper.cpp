#include <iostream>

#include "LustrineWrapper.hpp"

int main(void) {
	std::cout << "Wrapper experiment" << std::endl;
	Lustrine::SimulationParameters parameters;
	parameters.X = 30;
	parameters.Y = 20;
	parameters.Z = 30;

	Lustrine::Wrapper::Grid* grids = nullptr;
	Lustrine::Wrapper::allocate_grid_array(&grids, 2);

	Lustrine::Wrapper::Position position;

	Lustrine::Wrapper::allocate_grid(&grids[0], 10, 10, 10, true);
	Lustrine::Wrapper::init_grid_box(&parameters, &grids[1], 10, 10, 10, 0, {1.0, 0.0, 0.0, 1.0});

	Lustrine::Wrapper::init_simulation(&parameters, grids, &position, 1);

	std::cout << "simulate for 1000 steps !" << std::endl;

	float dt = 0.0;
	for (int i = 0; i < 1000; i++) {
		Lustrine::Wrapper::simulate(dt);
		dt += 0.01;
	}

	std::cout << "dt = " << dt << std::endl;

	Lustrine::Wrapper::free_grid(grids);
	Lustrine::Wrapper::cleanup_simulation();

}