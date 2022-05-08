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


void init_grid_from_magika_voxel(Grid* grid, const std::string& path, glm::vec3 position, MaterialType type) {
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
    grid->cells = std::vector<int>(grid->num_grid_cells, 0);
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
                grid->cells[grid_index] = color_index;
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
    grid->dynamic_solid = false;

    grid->position = position;
}

JsonWriter &operator<<(JsonWriter &j, const ogt_vox_transform &t) {
    // ogt_vox_transform matrix is in col-major order, so we transpose it here
    j.bArray();
    j << t.m00 << t.m10 << t.m20 << t.m30;
    j << t.m01 << t.m11 << t.m21 << t.m31;
    j << t.m02 << t.m12 << t.m22 << t.m32;
    j << t.m03 << t.m13 << t.m23 << t.m33;
    j.eArray();
    return j;
}

JsonWriter &operator<<(JsonWriter &j, const ogt_vox_rgba &c) {
    j.bObject();
    j.safeKey("r") << c.r;
    j.safeKey("g") << c.g;
    j.safeKey("b") << c.b;
    j.safeKey("a") << c.a;
    j.eObject();
    return j;
}

std::string read_vox_scene_json(const uint8_t *buffer, int64_t size) {
    std::stringstream s;
    const ogt_vox_scene *scene = ogt_vox_read_scene(buffer, size);
    if (!scene) throw std::runtime_error("failed to read scene");
    JsonWriter j{s};
    j.bObject();
    j.safeKey("models").bArray();
    int maxMat = 0;
    for (int m = 0; m < scene->num_models; ++m) {
        const ogt_vox_model *model = scene->models[m];
        j.bObject();
        j.safeKey("id") << m;
        j.safeKey("sizeX") << model->size_x;
        j.safeKey("sizeY") << model->size_y;
        j.safeKey("sizeZ") << model->size_z;
        j.safeKey("data").bArray();
        // color indices in x -> y -> z order. a color index of 0 means empty, all other indices mean solid and can be
        // used to index the scene's palette to obtain the color for the voxel
        uint32_t model_size = model->size_x * model->size_y * model->size_z;
        for (int i = 0; i < model_size; ++i) {
            j << model->voxel_data[i];
            maxMat = std::max<int>(maxMat, model->voxel_data[i]);
        }
        j.eArray();
        j.eObject();
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
        j.safeKey("parentGroup") << group.parent_group_index;
        j.safeKey("layer") << group.layer_index;
        j.safeKey("hidden") << group.hidden;
        j.eObject();
    }
    j.eArray();
    j.safeKey("palette").bArray();
    for (int p = 0; p <= maxMat; ++p) {
        j << scene->palette.color[p];
    }
    j.eArray();
    j.safeKey("materials").bArray();
    for (int m = 0; m <= maxMat; ++m) {
        const auto &mat = scene->materials.matl[m];
        j.bObject();
        j.safeKey("contentFlags") << mat.content_flags;
        j.safeKey("type") << mat.type;
        if (mat.content_flags & k_ogt_vox_matl_have_metal) j.safeKey("metal") << mat.metal;
        if (mat.content_flags & k_ogt_vox_matl_have_rough) j.safeKey("rough") << mat.rough;
        if (mat.content_flags & k_ogt_vox_matl_have_spec) j.safeKey("spec") << mat.spec;
        if (mat.content_flags & k_ogt_vox_matl_have_ior) j.safeKey("ior") << mat.ior;
        if (mat.content_flags & k_ogt_vox_matl_have_att) j.safeKey("att") << mat.att;
        if (mat.content_flags & k_ogt_vox_matl_have_flux) j.safeKey("flux") << mat.flux;
        if (mat.content_flags & k_ogt_vox_matl_have_emit) j.safeKey("emit") << mat.emit;
        if (mat.content_flags & k_ogt_vox_matl_have_ldr) j.safeKey("ldr") << mat.ldr;
        if (mat.content_flags & k_ogt_vox_matl_have_trans) j.safeKey("trans") << mat.trans;
        if (mat.content_flags & k_ogt_vox_matl_have_alpha) j.safeKey("alpha") << mat.alpha;
        if (mat.content_flags & k_ogt_vox_matl_have_d) j.safeKey("d") << mat.d;
        if (mat.content_flags & k_ogt_vox_matl_have_sp) j.safeKey("sp") << mat.sp;
        if (mat.content_flags & k_ogt_vox_matl_have_g) j.safeKey("g") << mat.g;
        if (mat.content_flags & k_ogt_vox_matl_have_media) j.safeKey("media") << mat.media;
        j.eObject();
    }
    j.eArray();
    j.eObject();
    ogt_vox_destroy_scene(scene);
    return s.str();
}
}