#pragma once

#include "Simulation.hpp"

namespace Lustrine {

	extern void simulate_fluid(Simulation* simulation, float dt);
	extern void simulate_sand(Simulation* simulation, float dt);
	extern void simulate_sand_v1(Simulation* simulation, float dt);

}