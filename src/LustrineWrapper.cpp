#include "LustrineWrapper.hpp"

#include <iostream>

namespace Lustrine {
	Simulation* simulation;
	void init_simulation_wrapper(const SimulationParameters* parameters) {
		std::cout << "wrapper init called" << std::endl;
	}
}