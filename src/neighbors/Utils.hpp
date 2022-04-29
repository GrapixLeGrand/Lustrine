#pragma once

#include "Simulation.hpp"

namespace Lustrine {

glm::vec3 get_cell_id_comp(const Simulation* simulation, glm::vec3 position, int i);
int get_cell_id(const Simulation* simulation, glm::vec3 position);

}