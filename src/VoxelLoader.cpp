#include <iostream>
#include <sstream>

#include "VoxelLoader.hpp"
#include "Lustrine.hpp"

#define OGT_VOX_IMPLEMENTATION
#include "../thirdparty/opengametools/ogt_vox.h"
#include "JsonWriter.h"


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


void init_grid_from_magika_voxel(Grid* grid, const std::string& path, MaterialType type) {
    const ogt_vox_scene* scene = load_vox_scene(path.c_str());

    if (scene->num_models < 1) {
        std::cout << "Emtpy voxel model" << std::endl;
    }
    std::cout << scene->num_layers << std::endl;
    const ogt_vox_model* model = scene->models[0];

    grid->X = model->size_x;
    grid->Y = model->size_y;
    grid->Z = model->size_z;

    grid->type = type;
    grid->num_grid_cells = model->size_x * model->size_y * model->size_z;
    grid->cells = std::vector<bool>(grid->num_grid_cells, false);
    grid->colors = std::vector<glm::vec4>(grid->num_grid_cells, {0, 0, 0, 0});
    grid->has_one_color_per_cell = true;

    int counter = 0;
    for (int x = 0; x < grid->X; x++) {
        for (int y = 0; y < grid->Y; y++) {
            for (int z = 0; z < grid->Z; z++) {
                int voxel_index = x + (y * model->size_x) + (z * model->size_x * model->size_y);
                int grid_index = (x * grid->Y * grid->Z) + (y * grid->Z) + z;
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
                grid_color /= 255.0;

                counter++;
            }
        }
    }

    grid->num_occupied_grid_cells = counter;
    ogt_vox_destroy_scene(scene);

    grid->sparse_solid = true;
    grid->dynamic_solid = true;

}

JsonWriter &operator<<(JsonWriter &j, const ogt_vox_transform &t) {
    j.bArray();
    j << t.m00 << t.m01 << t.m02 << t.m03;
    j << t.m10 << t.m11 << t.m12 << t.m13;
    j << t.m20 << t.m21 << t.m22 << t.m23;
    j << t.m30 << t.m31 << t.m32 << t.m33;
    j.eArray();
    return j;
}

std::string read_vox_scene_json(const uint8_t *buffer, int64_t size) {
    std::stringstream s;
    const ogt_vox_scene *scene = ogt_vox_read_scene(buffer, size);
    JsonWriter j{s};
    j.bObject();
    j.safeKey("models").bArray();
    for (int m = 0; m < scene->num_models; ++m) {
        const ogt_vox_model *model = scene->models[m];
        j.safeKey("id") << m;
        j.safeKey("size_x") << model->size_x;
        j.safeKey("size_y") << model->size_y;
        j.safeKey("size_z") << model->size_z;
        j.safeKey("data").bArray();
        uint32_t model_size = model->size_x * model->size_y * model->size_z;
        for (int i = 0; i < model_size; ++i) j << model->voxel_data[i];
        j.eArray();
    }
    j.eArray();
    j.safeKey("instances").bArray();
    for (int i = 0; i < scene->num_instances; ++i) {
        const ogt_vox_instance &inst = scene->instances[i];
        j.bObject();
        j.safeKey("id") << i;
        j.safeKey("name") << inst.name;
        j.safeKey("transform") << inst.transform;
        j.safeKey("model") << inst.model_index;
        j.safeKey("layer") << inst.layer_index;
        j.safeKey("group") << inst.group_index;
        j.safeKey("hidden") << inst.hidden;
        j.eObject();
    }
    j.eArray();
    j.safeKey("layers").bArray();
    for (int l = 0; l < scene->num_layers; ++l) {
        const ogt_vox_layer &layer = scene->layers[l];
        j.bObject();
        j.safeKey("id") << l;
        j.safeKey("name") << layer.name;
        j.safeKey("hidden") << layer.hidden;
        j.eObject();
    }
    j.eArray();
    j.safeKey("groups").bArray();
    for (int g = 0; g < scene->num_groups; ++g) {
        const ogt_vox_group &group = scene->groups[g];
        j.bObject();
        j.safeKey("id") << g;
        j.safeKey("transform") << group.transform;
        j.safeKey("parent_group") << group.parent_group_index;
        j.safeKey("layer") << group.layer_index;
        j.safeKey("hidden") << group.hidden;
        j.eObject();
    }
    j.eArray();
    j.eObject();
    ogt_vox_destroy_scene(scene);
    return s.str();
}
}