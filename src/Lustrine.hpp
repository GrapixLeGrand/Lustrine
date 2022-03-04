#pragma once

#include "Simulation.hpp"
#include "Kernels.hpp"
#include "Neighbors.hpp"
#include <vector>
#include "glm/glm.hpp"

//https://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf

namespace Lustrine {
struct Simulation;


extern void init_sim(Simulation* simulation, int particlesX, int particlesY, int particlesZ);
extern void fill_grid(Simulation* simulation);

extern float s_coor(const Simulation* simulation, float rl);

extern void simulate(Simulation* simulation, float dt);

};