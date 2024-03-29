#include "Utils.hpp"

namespace Lustrine {
    
    inline glm::vec3 get_cell_id_comp(const Simulation* simulation, glm::vec3 position, int i) {
        //glm::vec3 position = simulation->positions_star[i];
        position = glm::clamp(position, glm::vec3(simulation->cell_size * 0.5), glm::vec3(simulation->domainX - simulation->cell_size * 0.5, simulation->domainY - simulation->cell_size * 0.5, simulation->domainZ - simulation->cell_size * 0.5));
        position /= simulation->cell_size;
        position.x = (int) position.x;
        position.y = (int) position.y;
        position.z = (int) position.z;
        return position;
    }
    
}