#include <iostream>

#include "LustrineWrapper.hpp"

int main(void) {
	std::cout << "Wrapper experiment" << std::endl;
	Lustrine::SimulationParameters parameters;
	Lustrine::Wrapper::Grid grid;
	Lustrine::Wrapper::Position position;

	Lustrine::Wrapper::init_simulation_wrapper(&parameters, &grid, &position);
}