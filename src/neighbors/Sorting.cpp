
#include "Sorting.hpp"
#include "Utils.hpp"
#include "Simulation.hpp"
#include <string.h>

namespace Lustrine {
    namespace Sorting {
          
        void counting_sort(
            int* counts,
            int* initial_indices,
            const size_t num_positions,
            const size_t num_cells,
            int* sorted_indices
        ) {

            memset(counts, 0, (num_cells + 1) * sizeof(int));
            for (int i = 0; i < num_positions; i++) {
                counts[initial_indices[i]]++;
            }

            for (int i = 1; i <= num_cells; i++) {
                counts[i] += counts[i - 1];
            }

            for (int i = num_positions - 1; i >= 0; i--) {
                sorted_indices[counts[initial_indices[i]] - 1] = i;
                counts[initial_indices[i]]--;
            }

        }

    }
}