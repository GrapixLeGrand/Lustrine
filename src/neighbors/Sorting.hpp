#pragma once

#include "Lustrine_Export.h"
#include <cstddef>

namespace Lustrine {
    namespace Sorting {
        /**
        * @param counts: [0, num_cells - 1]
        * @param initial_indices: [0, num_positions - 1]
        * @param sorted_indices: [0, num_positions - 1]
        */
        extern void LUSTRINE_EXPORT counting_sort(
            int* counts,
            int* initial_indices,
            const size_t num_positions,
            const size_t num_cells,
            int* sorted_indices
        );

        extern void LUSTRINE_EXPORT counting_sort2(
            int* counts,
            int* initial_indices,
            const size_t num_positions,
            const size_t num_cells,
            int* sorted_indices
        );
    }
}