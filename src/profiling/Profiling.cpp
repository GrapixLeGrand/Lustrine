#include "Profiling.hpp"
#include <cassert>

namespace Lustrine {
	namespace Profiling {

		long long raw_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];
		double raw_durations[LUSTRINE_MAX_NUM_MEASUREMENTS];

		long get_num_observation() {
			return LUSTRINE_MAX_NUM_MEASUREMENTS;
		}

		long long get_cycles(int index) {
			assert(index >= 0 & index < LUSTRINE_MAX_NUM_MEASUREMENTS);
			return raw_cycles[index];
		}

		double get_duration(int index) {
			assert(index >= 0 & index < LUSTRINE_MAX_NUM_MEASUREMENTS);
			return raw_durations[index];
		}

	}
}