#pragma once

#include <iostream>
#include "Lustrine_Export.h"

#ifdef PLATFORM_WINDOWS
#define NOMINMAX
#include "Windows.h"
#elif PLATFORM_UNIX
#include <string.h>
#include "tsc_x86.h"
#endif

#define LUSTRINE_MAX_NUM_MEASUREMENTS 7

namespace Lustrine {
	
	namespace Profiling {
		//extern and cannot be static because of the header
		extern long long raw_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];
		extern double raw_durations[LUSTRINE_MAX_NUM_MEASUREMENTS];

#ifdef PLATFORM_WINDOWS

		static LARGE_INTEGER frequency;
		static LARGE_INTEGER start_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];
		
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
			//std::cout << (long long)start_cycles[index].QuadPart << std::endl;
		}

		static inline void stop_counter(int index) {
			LARGE_INTEGER tmp;
			QueryPerformanceCounter(&tmp);
			LONGLONG sub = tmp.QuadPart - start_cycles[index].QuadPart;
			raw_cycles[index] = (long long) sub;
			//std::cout << raw_cycles[index] << std::endl;
			raw_durations[index] = ((double)sub / frequency.QuadPart);
		}
#elif PLATFORM_UNIX
		#define QUENTIN_FREQUENCY 2700000000.0 //hz
		static unsigned long long start_cycles[LUSTRINE_MAX_NUM_MEASUREMENTS];

		static void init_profiling() {
			init_tsc();
			memset(start_cycles, 0, LUSTRINE_MAX_NUM_MEASUREMENTS * sizeof(unsigned long long));
			memset(raw_cycles, 0, LUSTRINE_MAX_NUM_MEASUREMENTS * sizeof(long long));
			memset(raw_durations, 0, LUSTRINE_MAX_NUM_MEASUREMENTS * sizeof(double));
		}

		static inline void start_counter(int index) {
			start_cycles[index] = start_tsc();
		}

		static inline void stop_counter(int index) {
			myInt64 result = stop_tsc(start_cycles[index]);
			raw_cycles[index] = result;
			raw_durations[index] = result / QUENTIN_FREQUENCY;
		}
#else
		static void init_profiling() { std::cout << "Lustrine::Profiling: disabled" << std::endl; }
		static inline void start_counter(int index) {}
		static inline void stop_counter(int index) {}
#endif
		extern "C" LUSTRINE_EXPORT long get_num_observation();
		extern "C" LUSTRINE_EXPORT long long get_cycles(int index);
		extern "C" LUSTRINE_EXPORT double get_duration(int index);

	}

	
}