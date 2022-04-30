#include <iostream>

#include "LustrineWrapper.hpp"

int main(void) {
	std::cout << "Wrapper experiment" << std::endl;
	Lustrine::SimulationParameters parameters;
	parameters.X = 30;
	parameters.Y = 20;
	parameters.Z = 30;

	Lustrine::Wrapper::SimulationData data;
	Lustrine::Wrapper::Grid sand_grids[1];
	Lustrine::Wrapper::Vec3 sand_positions[1];

	Lustrine::Wrapper::Grid solid_grids[2];
	Lustrine::Wrapper::Vec3 solid_positions[2];

	sand_positions[0] = {0, 0, 0};
	solid_positions[0] = { 10.0f, 0, 10.0f };
	solid_positions[1] = { 0, 0, 0 };

	Lustrine::Wrapper::init_grid_box(&parameters, &sand_grids[0], 10, 10, 10, { 0, 0, 0 }, { 1.0f, 0.0f, 0.0f, 1.0f }, 1);
	Lustrine::Wrapper::init_grid_box(&parameters, &solid_grids[0], 1, 5, 1, { 10.0f, 0, 10.0f }, { 1.0f, 0.0f, 0.0f, 1.0f }, 0);
	Lustrine::Wrapper::init_grid_box(&parameters, &solid_grids[1], 30, 1, 30, { 0, 0, 0 }, { 1.0f, 0.0f, 0.0f, 1.0f }, 0);

	Lustrine::Wrapper::init_simulation(
		&parameters,
		&data,
		sand_grids,
		1,
		solid_grids,
		2,
		1
	);

	int box = Lustrine::Wrapper::add_box({5.0f, 2.0f, 5.0f}, true, {0.5f, 0.5f, 0.5f});
	//int box_index = Lustrine::Bullet::add_box(bulletPhysic, {15, 15, 15}, true);

	int num_bodies = Lustrine::Wrapper::get_num_bodies();
	printf("Num bodies registered %d\n", num_bodies);
	std::cout << "simulate for 1000 steps !" << std::endl;

	int* indices = new int[num_bodies];
	memset(indices, -1, num_bodies);
	int actual = -1;
	Lustrine::Wrapper::check_collisions(box, indices, &actual);
	bool collide2 = Lustrine::Wrapper::do_collide(box);

	float dt = 0.01;
	Lustrine::Wrapper::Vec3 v {0, 0, 0};
	for (int i = 0; i < 1000; i++) {
		//std::cout << "simulate ! "<< std::endl;
		Lustrine::Wrapper::Vec3 v = Lustrine::Wrapper::get_position(box);
		//Lustrine::Wrapper::set_velocity(box, {0, 1, 0});
		printf("{%f, %f, %f}\n", v.x, v.y, v.z);
		Lustrine::Wrapper::simulate(dt, v, false, false);
		Lustrine::Wrapper::check_collisions(box, indices, &actual);

		if (actual > 0) {
			int a[5] = {indices[0], indices[1], indices[2], indices[3], indices[4] };
			bool collide = Lustrine::Wrapper::do_collide(box);
			printf("youhou\n");
		}

		//std::cout << particles_positions[33] << " " << particles_positions[34] << " " << particles_positions[35] << std::endl;
	}

	std::cout << "dt = " << dt << std::endl;
	Lustrine::Wrapper::cleanup_simulation();

}