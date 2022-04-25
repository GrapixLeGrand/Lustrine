#include "Profiling.hpp"

namespace Lustrine {

	long get_num_observation() {
		return LUSTRINE_MAX_NUM_MEASUREMENTS;
	}
	/*
	void get_cycles_observation(int index) {
		return Profiling::raw_cycles[index];
	}*/
}