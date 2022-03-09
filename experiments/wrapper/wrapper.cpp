#include <iostream>

#include "LustrineWrapper.hpp"

int main(void) {
	std::cout << "Wrapper experiment" << std::endl;
	Lustrine::SimulationParameters parameters;
	parameters.X = 30;
	parameters.Y = 20;
	parameters.Z = 30;

	Lustrine::Wrapper::Grid grid;
	Lustrine::Wrapper::Position position;

	Lustrine::Wrapper::allocate_grid(&grid, 10, 10, 10, true);
	Lustrine::Wrapper::init_simulation_wrapper(&parameters, &grid, &position, 1);
	Lustrine::Wrapper::free_grid(&grid);
	Lustrine::Wrapper::cleanup_simulation();

}