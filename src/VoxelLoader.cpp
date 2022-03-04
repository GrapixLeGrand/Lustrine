
#include "VoxelLoader.hpp"
#include <iostream>
#include "../thirdparty/opengametools/ogt_vox.h"
#define OGT_VOX_IMPLEMENTATION

namespace Lustrine {
    // a helper function to load a magica voxel scene given a filename.
const ogt_vox_scene* load_vox_scene(const char* filename, uint32_t scene_read_flags = 0)
{
    // open the file
#if defined(_MSC_VER) && _MSC_VER >= 1400
    FILE * fp;
    if (0 != fopen_s(&fp, filename, "rb"))
        fp = 0;
#else
    FILE * fp = fopen(filename, "rb");
#endif
    if (!fp)
        return NULL;

    // get the buffer size which matches the size of the file
    fseek(fp, 0, SEEK_END);
    uint32_t buffer_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    // load the file into a memory buffer
    uint8_t * buffer = new uint8_t[buffer_size];
    fread(buffer, buffer_size, 1, fp);
    fclose(fp);

    // construct the scene from the buffer
    const ogt_vox_scene * scene = ogt_vox_read_scene_with_flags(buffer, buffer_size, scene_read_flags);

    // the buffer can be safely deleted once the scene is instantiated.
    delete[] buffer;

    return scene;
}


void init_grid_from_magika_voxel(Grid* grid, const std::string& path) {
    const ogt_vox_scene* scene = load_vox_scene(path.c_str());

    if (scene->num_models < 1) {
        std::cout << "Emtpy voxel model" << std::endl;
    }

    const ogt_vox_model* model = scene->models[0];

    grid->X = model->size_x;
    grid->Y = model->size_y;
    grid->Z = model->size_z;

    grid->num_grid_cells = model->size_x * model->size_y * model->size_z;
    grid->cells = std::vector<bool> (grid->num_grid_cells, 0);
    grid->colors = std::vector<glm::vec4> (grid->num_grid_cells, {0, 0, 0, 0});
    grid->has_one_color_per_cell = true;

    int counter = 0;
    for (int x = 0; x < grid->X; x++) {
        for (int y = 0; y < grid->Y; y++) {
            for (int z = 0; z < grid->Z; z++) {
                int voxel_index = x + (y * model->size_x) + (z * model->size_x * model->size_y);
                int grid_index = (y * grid->X * grid->Z) + (x * grid->Z) + z;
                uint8_t color_index = model->voxel_data[voxel_index];

                if (color_index == 0) { //voxel is non existent
                    continue;
                }

                ogt_vox_rgba voxel_color = scene->palette.color[color_index];
                grid->cells[grid_index] = true;
                glm::vec4& grid_color = grid->colors[grid_index];
                grid_color.r = voxel_color.r;
                grid_color.g = voxel_color.g;
                grid_color.b = voxel_color.b;
                grid_color.a = voxel_color.a;

                counter++;
            }
        }
    }

    grid->num_occupied_grid_cells = counter;
    ogt_vox_destroy_scene(scene);
}

}