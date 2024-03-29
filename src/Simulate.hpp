#pragma once

#include "Simulation.hpp"
#include "Lustrine_Export.h"

namespace Lustrine {

	/**
	 * currently only simulate_sand and simulate_sand_v3 are stable
	 */

	extern LUSTRINE_EXPORT void simulate_fluid(Simulation* simulation, float dt);
	extern LUSTRINE_EXPORT void simulate_sand(Simulation* simulation, float dt);
	extern LUSTRINE_EXPORT void simulate_sand_credits(Simulation* simulation, float dt);
	extern LUSTRINE_EXPORT void simulate_sand_v1(Simulation* simulation, float dt);
	extern LUSTRINE_EXPORT void simulate_sand_v2(Simulation* simulation, float dt);
	extern LUSTRINE_EXPORT void simulate_sand_v3(Simulation* simulation, float dt);
	
}