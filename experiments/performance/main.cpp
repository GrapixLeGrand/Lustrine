#include <iostream>
//#include "tsc_x86.h"
//#include <profileapi.h>
#include <Windows.h>

/** Use to init the clock */
#define TIMER_INIT \
    LARGE_INTEGER frequency; \
    LARGE_INTEGER t1,t2; \
    double elapsedTime; \
    QueryPerformanceFrequency(&frequency);


/** Use to start the performance timer */
#define TIMER_START QueryPerformanceCounter(&t1);

/** Use to stop the performance timer and output the result to the standard stream. Less verbose than \c TIMER_STOP_VERBOSE */
#define TIMER_STOP \
    QueryPerformanceCounter(&t2); \
    elapsedTime=(float)(t2.QuadPart-t1.QuadPart)/frequency.QuadPart; \
    std::cout<< elapsedTime <<std::endl;


int main(int argc, char** args) {
	
	std::cout << "performance observed" << std::endl;
	TIMER_INIT;
	int lim;
	std::cin >> lim;
	TIMER_START;
	float a = 0;
	for (int i = 0; i < lim; i++) {
		a += i;
	}
	TIMER_STOP;

	std::cout << a << std::endl;
	return 0;
}