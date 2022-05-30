#pragma once

#include "Simulation.hpp"

namespace Lustrine {

struct CellIndices {
    int x;
    int y;
    int z;
};

static inline CellIndices get_cell_indices(const Simulation* simulation, glm::vec3 position) {
    position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
    position /= simulation->cell_size;
    CellIndices result;
    result.x = (int) position.x;
    result.y = (int) position.y;
    result.z = (int) position.z;
    return result;
}

inline glm::vec3 get_cell_id_comp(const Simulation* simulation, glm::vec3 position, int i);
inline static int get_cell_id(const Simulation* simulation, glm::vec3 position) {
    //position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5f), glm::vec3(simulation->domainX - simulation->cell_size * 0.5f, simulation->domainY - simulation->cell_size * 0.5f, simulation->domainZ - simulation->cell_size * 0.5f));
    position /= simulation->cell_size;
    int cell_id =
            ((int) position.y) * simulation->gridX * simulation->gridZ +
            ((int) position.x) * simulation->gridZ +
            ((int) position.z); 
    //cell_id = glm::clamp(cell_id, 0, simulation->num_grid_cells - 1);
    return cell_id;
}


}