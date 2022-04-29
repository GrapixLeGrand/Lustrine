#include "Profiling.hpp"
#include <cassert>

//bellow is to allow tsc 86 header to be included in other header files
#ifdef PLATFORM_UNIX
void init_tsc() {
	; // no need to initialize anything for x86
}

myInt64 start_tsc(void) {
    tsc_counter start;
    CPUID();
    RDTSC(start);
    return COUNTER_VAL(start);
}

myInt64 stop_tsc(myInt64 start) {
	tsc_counter end;
	RDTSC(end);
	CPUID();
	return COUNTER_VAL(end) - start;
}
#endif

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