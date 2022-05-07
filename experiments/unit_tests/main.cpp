#include <iostream>
#include <vector>
#include "Lustrine.hpp"
#include "profiling/Profiling.hpp"
#include "Simulate.hpp"
#include <string>
#include "neighbors/Sorting.hpp"
#include "Lustrine.hpp"
#include <cstdlib>
#include <cassert>
#include <algorithm>

std::vector<int> generate_random_ints(size_t num, int max) {
	std::vector<int> result(num);

	for (size_t i = 0; i < num; i++) {
		result[i] = rand() % max;
	}
	return result;
}

void print_int_array(int* arr, size_t size) {
	for (size_t i = 0; i < size; i++) {
		std::cout << arr[i] << " ";
	}
	std::cout << std::endl;
}

void permute(int* permute, size_t size, int* permutation_indices) {

	int* tmp = (int*) malloc(size * sizeof(int));
	memcpy(tmp, permute, size * sizeof(int));

	print_int_array(tmp, size);
	print_int_array(permute, size);

	for (int i = 0; i < size; i++) {
		int permuted_index = permutation_indices[i];
		int real_index = tmp[permuted_index];

		permute[i] = tmp[permutation_indices[i]];
	}
	free(tmp);
}

void test_1() {

	std::cout << "Test 1 start" << std::endl;

	size_t num = 20;
	size_t num_cells = 5;
	std::vector<int> cell_indices = generate_random_ints(num, num_cells);

	int* cell_indices_to_sort = (int*)malloc(num * sizeof(int));
	std::vector<int> counts(num_cells);

	std::cout << "initial: ";
	print_int_array(cell_indices.data(), num);

	int* tmp = (int*)malloc(num * sizeof(int));
	memcpy(tmp, cell_indices.data(), num * sizeof(int));

	int* counts_ptr = counts.data();
	int* initial = cell_indices.data();
	Lustrine::Sorting::counting_sort(counts_ptr, initial, num, num_cells, cell_indices_to_sort);
	std::sort(cell_indices.begin(), cell_indices.end());

	//std::cout << "permutation: ";
	//print_int_array(cell_indices_to_sort, num);

	//tmp contains the original input and must become sorted with the counting sort indices result
	permute(tmp, num, cell_indices_to_sort);

	//std::cout << "original sorted: ";
	//print_int_array(cell_indices.data(), num);
	//std::cout << "original permuted: ";
	//print_int_array(tmp, num);

	int result = memcmp(cell_indices.data(), tmp, num * sizeof(int));
	if (result == 0) {
		std::cout << "--> Success!" << "\n";
	}
	else {
		std::cout << "--> failure!" << "\n";
	}
	assert(result == 0);

	std::cout << "Test 1 end" << std::endl;
}


int main(int argc, char** args) {


	std::cout << "Unit testing" << std::endl;
	test_1();

	return 0;
}