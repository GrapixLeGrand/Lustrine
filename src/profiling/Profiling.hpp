#pragma once

#include <iostream>
#include "Lustrine_Export.h"

#ifdef PLATFORM_WINDOWS
#define NOMINMAX
#include "Windows.h"
#endif

#define LUSTRINE_MAX_NUM_MEASUREMENTS 3

namespace Lustrine {
	
	namespace Profiling {
		
#ifdef PLATFORM_WINDOWS

		static LARGE_INTEGER frequency;
		static LARGE_INTEGER start_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];
		static unsigned long long raw_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];
		static double raw_durations[LUSTRINE_MAX_NUM_MEASUREMENTS];

		static void init_profiling() {
			std::cout << "Lustrine::Profiling: enabled" << std::endl;
			QueryPerformanceFrequency(&frequency);
			for (size_t i = 0; i < LUSTRINE_MAX_NUM_MEASUREMENTS; i++) {
				raw_cycles[i] = 0;
				raw_durations[i] = 0.0;
			}
		}

		static inline void start_counter(int index) {
			QueryPerformanceCounter(&start_cycles[index]);
		}

		static inline void stop_counter(int index) {
			LARGE_INTEGER tmp;
			QueryPerformanceCounter(&tmp);
			LONGLONG sub = tmp.QuadPart - start_cycles[index].QuadPart;
			raw_cycles[index] = (unsigned long long) sub;
			raw_durations[index] = (double)sub / frequency.QuadPart;
		}

		extern "C" LUSTRINE_EXPORT long get_num_observation();
		extern "C" LUSTRINE_EXPORT void get_cycles_observation(int index);
#else
		static void init_profiling() { std::cout << "Lustrine::Profiling: disabled" << std::endl; }
		static inline void start_counter(int index) {}
		static inline void stop_counter(int index) {}
#endif
	}

	
}