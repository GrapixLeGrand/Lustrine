#include <iostream>
//#include "tsc_x86.h"
//#include <profileapi.h>
#include <iostream>
#include <vector>
#include "Lustrine.hpp"
#include "profiling/Profiling.hpp"
#include "Simulate.hpp"

void benchmark_single_iter(Lustrine::SimulationParameters* parameters, Lustrine::Simulation* simulation, int iter, Lustrine::Simulate_fun fun) {
	
	
	long long total = 0;
	std::vector<long long> cycles(Lustrine::Profiling::get_num_observation());
	std::vector<double> durations(Lustrine::Profiling::get_num_observation());

	std::vector<Lustrine::Grid> grids (1);
	Lustrine::init_grid_box(parameters, &grids[0], 20, 20, 20, {1.0, 1.0, 1.0}, {0, 0, 0, 1.0}, Lustrine::MaterialType::SAND);
	std::vector<Lustrine::Grid> solid_grids(1);
	Lustrine::init_grid_box(parameters, &solid_grids[0], 20, 20, 20, { 20, 1, 20 }, { 0, 0, 0, 1.0 }, Lustrine::MaterialType::SOLID);

	Lustrine::init_simulation(parameters, simulation, grids, solid_grids);

	simulation->simulate_fun = fun;

	for (int i = 0; i < iter; i++) {
		Lustrine::simulate(simulation, 0.01f);
		for (int j = 0; j < cycles.size(); j++) {
			cycles[j] += Lustrine::Profiling::get_cycles(j);
			durations[j] += Lustrine::Profiling::get_duration(j);
		}
	}

	for (int i = 0; i < cycles.size(); i++) {
		std::cout << i << ": " << (((double)cycles[i]) / iter) << " cylces\t" << (((double)durations[i]) / ((double)iter)) << " ms" << std::endl;
	}

	Lustrine::clean_simulation(simulation);
}

int main(int argc, char** args) {


	std::cout << "Performance profiling" << std::endl;


	Lustrine::Simulation simulation;
	Lustrine::SimulationParameters parameters;
	parameters.X = 50;
	parameters.Y = 40;
	parameters.Z = 50;


	std::cout << "Inititalized simulation" << std::endl;
	
	//warm up
	benchmark_single_iter(&parameters, &simulation, 10, Lustrine::simulate_sand);
	benchmark_single_iter(&parameters, &simulation, 10, Lustrine::simulate_sand_v1);

	return 0;
}