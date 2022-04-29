
#include "Sorting.hpp"
#include "Utils.hpp"
#include "Simulation.hpp"

namespace Lustrine {
    namespace Sorting {
        
        void counting_sort(Simulation* simulation) {

            std::fill(simulation->counts.begin(), simulation->counts.end(), 0);
            std::fill(simulation->particle_cell_index_to_index.begin(), simulation->particle_cell_index_to_index.end(), std::make_pair(0, 0));

            for (int i = 0; i < simulation->num_particles; i++) {
                int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
                simulation->counts[cell_id]++;
            }

            for (int i = 1; i < simulation->counts.size(); i++) {
                simulation->counts[i] += simulation->counts[i - 1];
            }

            for (int i = simulation->num_particles - 1; i >= 0; i--) {
                int cell_id = get_cell_id(simulation, simulation->positions_star[i]);
                simulation->counts[cell_id]--;
                simulation->particle_cell_index_to_index[simulation->counts[cell_id]].first = cell_id;
                simulation->particle_cell_index_to_index[simulation->counts[cell_id]].second = i;
                simulation->positions_star_copy[simulation->counts[cell_id]] = simulation->positions_star[i];
            }

        }

    }
}