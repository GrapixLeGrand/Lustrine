#include "LustrineWrapper.hpp"

#include <iostream>

namespace Lustrine {
	Simulation* simulation;
	void init_simulation_wrapper() {
		std::cout << "wrapper init called" << std::endl;
	}
}