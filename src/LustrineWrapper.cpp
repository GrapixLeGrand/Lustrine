#include "LustrineWrapper.hpp"

#include <iostream>

namespace Lustrine {
namespace Wrapper {

	//Simulation* simulation;
	void init_simulation_wrapper(const SimulationParameters* parameters, Grid* grids, Position* positions) {
		std::cout << "wrapper init called" << std::endl;
	}

	void init_grid_box_wrapper(const SimulationParameters* parameters, Grid* grid, int X, int Y, int Z, int type, Color color) {
		std::cout << "init called" << std::endl;
	}

	void simulate_wrapper(float dt) {
		std::cout << "dt is " << dt << std::endl;
	}

}
}