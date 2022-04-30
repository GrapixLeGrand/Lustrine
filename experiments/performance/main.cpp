#include <iostream>
//#include "tsc_x86.h"
//#include <profileapi.h>
#include <iostream>
#include <vector>
#include "Lustrine.hpp"
#include "profiling/Profiling.hpp"
#include "Simulate.hpp"
#include <string>
#include <fstream>

struct benchmark_result {
	std::vector<std::vector<double>> results;
	size_t num_measures = 0;
	size_t num_samples = 0;
};

const std::string out_file_name = "bench.csv";

benchmark_result benchmark_single_iter(Lustrine::SimulationParameters* parameters, Lustrine::Simulation* simulation, int iter, Lustrine::Simulate_fun fun) {
	
	
	long long total = 0;
	std::vector<long long> cycles(Lustrine::Profiling::get_num_observation());
	std::vector<double> durations(Lustrine::Profiling::get_num_observation());
	std::vector<std::vector<double>> results (iter, std::vector<double>(Lustrine::Profiling::get_num_observation(), 0.0));

	std::vector<Lustrine::Grid> grids (1);
	Lustrine::init_grid_box(parameters, &grids[0], 20, 20, 20, {1.0, 1.0, 1.0}, {0, 0, 0, 1.0}, Lustrine::MaterialType::SAND);
	std::vector<Lustrine::Grid> solid_grids(1);
	Lustrine::init_grid_box(parameters, &solid_grids[0], 20, 20, 20, { 20, 1, 20 }, { 0, 0, 0, 1.0 }, Lustrine::MaterialType::SOLID);

	Lustrine::init_simulation(parameters, simulation, grids, solid_grids);

	simulation->simulate_fun = fun;

	for (int i = 0; i < iter; i++) {
		Lustrine::simulate(simulation, 0.01f);
		for (int j = 0; j < cycles.size(); j++) {
			long long cycles_observed = Lustrine::Profiling::get_cycles(j);
			double duration_observed = Lustrine::Profiling::get_duration(j);
			results[i][j] = cycles_observed;
			cycles[j] += cycles_observed;
			durations[j] += duration_observed;
		}
	}

	for (int i = 0; i < cycles.size(); i++) {
		std::cout << std::fixed << i << ": " << (((double)cycles[i]) / iter) << " cylces\t" << (((double)durations[i]) / ((double)iter)) << " s" << std::endl;
	}

	Lustrine::clean_simulation(simulation);

	benchmark_result r = { results, (size_t) iter, (size_t) Lustrine::Profiling::get_num_observation() };

	return r;
}

void write_result_to_file(benchmark_result& result) {
	std::ifstream f(out_file_name.c_str());
    if (!f.good()) {
		std::cout << "create file" << std::endl;

	}
	
}

int main(int argc, char** args) {


	std::cout << "Performance profiling" << std::endl;


	Lustrine::Simulation simulation;
	Lustrine::SimulationParameters parameters;
	parameters.X = 50;
	parameters.Y = 40;
	parameters.Z = 50;
	parameters.particleRadius = 0.5f;
    parameters.particleDiameter = 1.0f;

	std::cout << "Inititalized simulation" << std::endl;
	int iter = 100;
	//warm up
	benchmark_single_iter(&parameters, &simulation, 1000, Lustrine::simulate_sand);

	benchmark_single_iter(&parameters, &simulation, iter, Lustrine::simulate_sand);
	benchmark_single_iter(&parameters, &simulation, iter, Lustrine::simulate_sand_v1);

	return 0;
}