#include <iostream>

#define OGT_VOX_IMPLEMENTATION
#include "opengametools/ogt_vox.h"
#include <cassert>
#include <algorithm>
#include <fstream>  


const ogt_vox_scene* load_vox_scene(const char* filename, uint32_t scene_read_flags = 0)
{
    // open the file
#if defined(_MSC_VER) && _MSC_VER >= 1400
    FILE* fp;
    if (0 != fopen_s(&fp, filename, "rb"))
        fp = 0;
#else
    FILE* fp = fopen(filename, "rb");
#endif
    if (!fp)
        return NULL;

    // get the buffer size which matches the size of the file
    fseek(fp, 0, SEEK_END);
    uint32_t buffer_size = ftell(fp);
    fseek(fp, 0, SEEK_SET);

    // load the file into a memory buffer
    uint8_t* buffer = new uint8_t[buffer_size];
    fread(buffer, buffer_size, 1, fp);
    fclose(fp);

    // construct the scene from the buffer
    const ogt_vox_scene* scene = ogt_vox_read_scene_with_flags(buffer, buffer_size, scene_read_flags);

    // the buffer can be safely deleted once the scene is instantiated.
    delete[] buffer;

    return scene;
}

void save_vox_scene(const char* pcFilename, const ogt_vox_scene* scene)
{
    // save the scene back out. 
    uint32_t buffersize = 0;
    uint8_t* buffer = ogt_vox_write_scene(scene, &buffersize);
    if (!buffer)
        return;

    // open the file for write
#if defined(_MSC_VER) && _MSC_VER >= 1400
    FILE* fp;
    if (0 != fopen_s(&fp, pcFilename, "wb"))
        fp = 0;
#else
    FILE* fp = fopen(pcFilename, "wb");
#endif
    if (!fp) {
        ogt_vox_free(buffer);
        return;
    }

    fwrite(buffer, buffersize, 1, fp);
    fclose(fp);
    ogt_vox_free(buffer);
}

#define NUM_CHAR_PATH 200
static char input_path[NUM_CHAR_PATH] = "C:/Users/miste/Desktop/gamelab2022-hyrule-team4/Source/game/Catrine/Content/vox_files/level2.vox";//LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_model.vox";
static char output_path[NUM_CHAR_PATH] = "C:/Users/miste/Desktop/gamelab2022-hyrule-team4/Source/game/Catrine/Content/vox_files/level2.vox";
double target_density = 0.01;
uint32_t subdivision = 1;

uint32_t count_voxels(const uint8_t* data, uint32_t size) {
    uint32_t counter = 0;
    for (uint32_t i = 0; i < size; i++) {
        if (data[i] > 0) counter++;
    }
    return counter;
}

uint32_t getsize(const ogt_vox_model* model) {
    return model->size_x * model->size_y * model->size_z;
}

int main(int argc, char** args) {
	
    if (argc > 1 && false) {
        std::cout << "using custom arguments" << "\n";
        assert(argc == 5);
        char* dummy = NULL;
        memcpy(input_path, args[1], strlen(args[1]));
        memcpy(output_path, args[2], strlen(args[2]));
        target_density = strtod(args[3], &dummy);
        subdivision = (uint32_t) strtol(args[4], &dummy, 10);
    }
    else {
        std::cout << "using default arguments" << "\n";
    }

    std::cout << "Arguments:" << "\n";
    std::cout << "\tin:\t" << input_path << "\n";
    std::cout << "\tout:\t" << output_path << "\n";
    std::cout << "\tdensity:\t" << target_density << "\n";
    std::cout << "\tsubdivision:\t" << subdivision << "\n";

	std::cout << "level maker start" << std::endl;
    std::cout << "model at path " << input_path << std::endl;//LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_model.vox" << std::endl;

	const ogt_vox_scene* scene = load_vox_scene(input_path);
    assert(scene != nullptr);

    std::cout << "scene properties:" << std::endl;
    std::cout << "\tnum models:\t" << scene->num_models << std::endl;
    std::cout << "\tnum groups:\t" << scene->num_groups << std::endl;
    std::cout << "\tnum instances:\t" << scene->num_instances << std::endl;
       
    assert(scene->num_models == 1);
    
    const ogt_vox_model* model = scene->models[0];
    
    std::cout << "Model properties:" << std::endl;
    std::cout << "\tx:\t" << model->size_x << std::endl;
    std::cout << "\ty:\t" << model->size_y << std::endl;
    std::cout << "\tz:\t" << model->size_z << std::endl;
    
    uint32_t voxel_per_subdivision = subdivision * subdivision * subdivision;

    uint32_t division_size_x = (model->size_x / subdivision) + (model->size_x % subdivision > 0 ? 1 : 0);
    uint32_t division_size_y = (model->size_y / subdivision) + (model->size_y % subdivision > 0 ? 1 : 0);
    uint32_t division_size_z = (model->size_z / subdivision) + (model->size_z % subdivision > 0 ? 1 : 0);
    
    std::cout << "\tsubdiv x:\t" << division_size_x << std::endl;
    std::cout << "\tsubdiv y:\t" << division_size_y << std::endl;
    std::cout << "\tsubdiv z:\t" << division_size_z << std::endl;
    
    uint32_t voxel_count = division_size_x * division_size_y * division_size_z;
    ogt_vox_model* subdivided_model = (ogt_vox_model*)_vox_calloc(sizeof(ogt_vox_model) + voxel_count);
    subdivided_model->size_x = division_size_x;
    subdivided_model->size_y = division_size_y;
    subdivided_model->size_z = division_size_z;
    subdivided_model->voxel_data = (uint8_t*)&subdivided_model[1]; // (uint8_t*)(((char*)subdivided_model) + sizeof(ogt_vox_model) - sizeof(uint8_t*));

    memset((void*) subdivided_model->voxel_data, 0, voxel_count);

    uint32_t occupied_voxel_count = 0;
    uint32_t total_counter = 0;
    for (uint32_t z = 0; z < division_size_z; z++) {
        for (uint32_t y = 0; y < division_size_y; y++) {
            for (uint32_t x = 0; x < division_size_x; x++) {

                uint32_t sub_max_x = std::min((x + 1) * subdivision, model->size_x);
                uint32_t sub_max_y = std::min((y + 1) * subdivision, model->size_y);
                uint32_t sub_max_z = std::min((z + 1) * subdivision, model->size_z);

                uint32_t current_voxel_volume = (sub_max_x % subdivision) * (sub_max_y % subdivision) * (sub_max_z % subdivision);

                uint32_t counter = 0;
                uint8_t last_palette_index = 0;
                   
                //std::cout << "z: " << sub_z << " " << sub_max_z << std::endl;

                for (uint32_t sub_z = z * subdivision; sub_z < sub_max_z; sub_z++) {
                    for (uint32_t sub_y = y * subdivision; sub_y < sub_max_y; sub_y++) {
                        for (uint32_t sub_x = x * subdivision; sub_x < sub_max_x; sub_x++) {
                        

                            uint32_t index = sub_x + (sub_y * model->size_x) + (sub_z * model->size_x * model->size_y);
                            if (model->voxel_data[index] != 0) {
                                counter++;
                                last_palette_index = model->voxel_data[index];
                            }
                            
                        }
                    }
                }
                
                total_counter += counter;
                double density = ((double)counter) / ((double)current_voxel_volume);
               
                if (density > target_density) {
                    occupied_voxel_count++;
                    uint32_t subdivided_index = x + (y * division_size_x) + (z * division_size_x * division_size_y);
                    ((uint8_t*)subdivided_model->voxel_data)[subdivided_index] = last_palette_index;
                }

            }
        }
    }

    uint32_t initial_dims = getsize(model);
    uint32_t initial_occupied = count_voxels(model->voxel_data, initial_dims);
    uint32_t initial_dims_sub = getsize(subdivided_model);
    uint32_t initial_occupied_sub = count_voxels(subdivided_model->voxel_data, initial_dims_sub);
    
    uint8_t* tmp_voxel_data = new uint8_t[subdivided_model->size_x * subdivided_model->size_y * subdivided_model->size_z];
    assert(tmp_voxel_data != nullptr);

    memcpy(tmp_voxel_data, subdivided_model->voxel_data, subdivided_model->size_x * subdivided_model->size_y * subdivided_model->size_z * sizeof(uint8_t));

    //here remove unseen voxels
    for (uint32_t z = 1; z < subdivided_model->size_z - 1; z++) {
        for (uint32_t y = 1; y < subdivided_model->size_y - 1; y++) {
            for (uint32_t x = 1; x < subdivided_model->size_x - 1; x++) {
               
                uint32_t X = subdivided_model->size_x;
                uint32_t Y = subdivided_model->size_y;

                uint32_t xx = x;
                uint32_t yy = y;
                uint32_t zz = z;

                if (
                    tmp_voxel_data[(x - 1) + y * X + z * X * Y] > 0 &&
                    tmp_voxel_data[(x + 1) + y * X + z * X * Y] > 0 &&
                    tmp_voxel_data[x + (y - 1) * X + z * X * Y] > 0 &&
                    tmp_voxel_data[x + (y + 1) * X + z * X * Y] > 0 &&
                    tmp_voxel_data[x + y * X + (z - 1) * X * Y] > 0 &&
                    tmp_voxel_data[x + y * X + (z + 1) * X * Y] > 0
                    ) {
                    ((uint8_t*)(subdivided_model->voxel_data))[x + y * X + z * X * Y] = 0; //remove the voxel if not visible
                }

            }
        }
    }

    delete[] tmp_voxel_data;
    subdivided_model->voxel_hash = _vox_hash(subdivided_model->voxel_data, subdivided_model->size_x * subdivided_model->size_y * subdivided_model->size_z);

    uint32_t initial_dims_simp = getsize(subdivided_model);
    uint32_t initial_occupied_simp = count_voxels(subdivided_model->voxel_data, initial_dims_sub);

    std::cout << "stats:" << "\n"
        << "\tinitial:\t" << "\ttotal vox\t" << initial_dims << "\tused:\t" << initial_occupied << "\toccupency:\t" << ((double)initial_occupied) / ((double)initial_dims) << "%\n"
        << "\tsubdiv:\t" << "\ttotal vox\t" << initial_dims_sub << "\tused:\t" << initial_occupied_sub << "\toccupency:\t" << ((double)initial_occupied_sub) / ((double)initial_dims_sub) << "%\n"
        << "\tsimplified:\t" << "\ttotal vox\t" << initial_dims_simp << "\tused:\t" << initial_occupied_simp << "\toccupency:\t" << ((double)initial_occupied_simp) / ((double)initial_dims_simp) << "%\n";

    //saving the subdivided model
    scene->models[0] = subdivided_model;
    save_vox_scene(output_path, scene);

    std::cout << "written scene" << std::endl;

    ogt_vox_free((void*)model);
    ogt_vox_destroy_scene(scene);
    std::cout << "level maker exit" << std::endl;
}