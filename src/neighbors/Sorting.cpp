
#include "Sorting.hpp"
#include "Utils.hpp"
#include "Simulation.hpp"
#include <string.h>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
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

        void counting_sort2(
            int* counts,
            int* initial_indices,
            const size_t num_positions,
            const size_t num_cells,
            int* sorted_indices
        ) {
            std::unordered_set<int> used_cells;
            //memset(counts, 0, (num_cells + 1) * sizeof(int));
            std::unordered_map<int, int> counts_map;
            for (int i = 0; i < num_positions; i++) {
                counts_map[initial_indices[i]]++;
                used_cells.insert(initial_indices[i]);
            }

            std::vector<int> used_cells_vec;
            used_cells_vec.reserve(used_cells.size());
            for (auto used_cell: used_cells) {
                used_cells_vec.push_back(used_cell);
            }
            std::sort(used_cells_vec.begin(), used_cells_vec.end());
            int current_count = 0;

            for (const int& used_cell : used_cells_vec) {
                current_count += counts_map[used_cell];
                counts_map[used_cell] = current_count;
            }

            for (int i = num_positions - 1; i >= 0; i--) {
                sorted_indices[counts_map[initial_indices[i]] - 1] = i;
                counts_map[initial_indices[i]]--;
            }

        }

    }
}