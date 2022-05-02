#pragma once

#include "Simulation.hpp"

namespace Lustrine {

inline glm::vec3 get_cell_id_comp(const Simulation* simulation, glm::vec3 position, int i);
inline static int get_cell_id(const Simulation* simulation, glm::vec3 position) {

    position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
    position /= simulation->cell_size;
    int cell_id =
            ((int) position.y) * simulation->gridX * simulation->gridZ +
            ((int) position.x) * simulation->gridZ +
            ((int) position.z);
            
    return cell_id;
}


}