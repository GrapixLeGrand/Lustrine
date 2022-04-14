#include <iostream>

#define OGT_VOX_IMPLEMENTATION
#include "opengametools/ogt_vox.h"
#define OGT_VOXEL_MESHIFY_IMPLEMENTATION
#include "opengametools/ogt_voxel_meshify.h"
#include <cassert>
#include <algorithm>
#include <fstream>  

#include "glm-0.9.9.8/glm/glm.hpp"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm-0.9.9.8/glm/gtx/hash.hpp"
#include <unordered_map>

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

/**
* From opengametool
*/
bool write_mesh_to_fbx(const char* output_filename, const ogt_mesh* mesh, const char* mesh_name)
{
#if defined(_MSC_VER) && _MSC_VER >= 1400  
    FILE* fout = NULL;
    fopen_s(&fout, output_filename, "wb");
#else
    FILE* fout = fopen(output_filename, "wb");
#endif
    if (!fout) {
        return false;
    }

    fprintf(fout,
        "; FBX 6.1.0 project file\n"
        "; ----------------------------------------------------\n"
        "\n"
        "FBXHeaderExtension:  {\n"
        "\tFBXHeaderVersion: 1003\n"
        "\tFBXVersion: 6100\n"
        "\tCreator: \"http://github.com/jpaver/opengametools vox2fbx\"\n"
        "}\n"
        "\n"
        "; Object definitions\n"
        "; ------------------------------------------------------------------\n"
        "\n"
        "Definitions:  {\n"
        "\tVersion: 100\n"
        "\tCount: 1\n"
        "\tObjectType: \"Model\" {\n"
        "\t\tCount: 1\n"
        "\t}\n"
        "}\n"
        "\n"
        "; Object properties\n"
        "; ------------------------------------------------------------------\n"
        "\n"
        "Objects:  {\n"
    );

    // start the model.
    fprintf(fout,
        "\tModel: \"%s\", \"Mesh\" {\n"
        "\t\tVersion: 232\n",
        mesh_name);

    // write the vertices
    {
        fprintf(fout,
            "\t\tVertices:");
        for (uint32_t i = 0; i < mesh->vertex_count; i++) {
            fprintf(fout, "%s%f,%f,%f", ((i > 0) ? "," : ""), mesh->vertices[i].pos.x, mesh->vertices[i].pos.y, mesh->vertices[i].pos.z);
        }
        fprintf(fout, "\n");
    }

    // write the vertex indices
    {
        fprintf(fout,
            "\t\tPolygonVertexIndex: ");
        for (uint32_t i = 0; i < mesh->index_count; i += 3) {
            uint32_t i0 = mesh->indices[i + 2];
            uint32_t i1 = mesh->indices[i + 1];
            uint32_t i2 = mesh->indices[i + 0];
            fprintf(fout, "%s%u,%u,-%u", ((i > 0) ? "," : ""), i0, i1, (i2 + 1));
        }
        fprintf(fout, "\n");
    }
    fprintf(fout,
        "\t\tGeometryVersion: 124\n"
    );

    // write the vertex normals layer element
    {
        fprintf(fout,
            "\t\tLayerElementNormal: 0 {\n"
            "\t\t\tVersion: 101\n"
            "\t\t\tName: \"\"\n"
            "\t\t\tMappingInformationType: \"ByVertice\"\n"
            "\t\t\tReferenceInformationType: \"Direct\"\n"
        );
        // write colors array
        fprintf(fout,
            "\t\t\tNormals: "
        );
        for (uint32_t i = 0; i < mesh->vertex_count; i++) {
            float x = mesh->vertices[i].normal.x;
            float y = mesh->vertices[i].normal.y;
            float z = mesh->vertices[i].normal.z;
            // palette color
            fprintf(fout, "%s%f,%f,%f", ((i > 0) ? "," : ""), x, y, z);
        }
        fprintf(fout, "\n");
        fprintf(fout,
            "\t\t}\n");
    }

    // write the vertex colors layer element
    {
        fprintf(fout,
            "\t\tLayerElementColor: 0 {\n"
            "\t\t\tVersion: 101\n"
            "\t\t\tName: \"colorSet1\"\n"
            "\t\t\tMappingInformationType: \"ByPolygonVertex\"\n"
            "\t\t\tReferenceInformationType: \"Direct\"\n"
        );
        // write colors array
        fprintf(fout,
            "\t\t\tColors: "
        );
        for (uint32_t i = 0; i < mesh->index_count; i++) {
            uint32_t index = mesh->indices[i];
            float r = (mesh->vertices[index].color.r / 255.0f);
            float g = (mesh->vertices[index].color.g / 255.0f);
            float b = (mesh->vertices[index].color.b / 255.0f);
            float a = 1.0f;
            // palette color
            fprintf(fout, "%s%f,%f,%f,%f", ((i > 0) ? "," : ""), r, g, b, a);
        }
        fprintf(fout, "\n");

        fprintf(fout,
            "\t\t}\n");
    }
    // write the layers
    fprintf(fout,
        "\t\tLayer: 0 {\n"
        "\t\t\tVersion: 100\n"
        "\t\t\tLayerElement: {\n"
        "\t\t\t\tType: \"LayerElementNormal\"\n"
        "\t\t\t\tTypedIndex: 0\n"
        "\t\t\t}\n"
        "\t\t\tLayerElement: {\n"
        "\t\t\t\tType: \"LayerElementColor\"\n"
        "\t\t\t\tTypedIndex: 0\n"
        "\t\t\t}\n"
        "\t\t}\n"
    );

    // write the tail of the model 
    fprintf(fout,
        "\t}\n"
        "}\n"
        "\n");

    // write the connections
    fprintf(fout,
        "; Object connections\n"
        "; ------------------------------------------------------------------\n"
        "Connections:  {\n"
        "\tConnect: \"OO\", \"%s\", \"Model::Scene\"\n"
        "}\n",
        mesh_name);

    fclose(fout);
    return true;
}


#define NUM_CHAR_PATH 200
static char input_path[NUM_CHAR_PATH] = LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_model.vox";
static char output_path[NUM_CHAR_PATH] = "out.obj";


int main(char** args, int argc) {
	
    /*
    if (argc > 1) {
        std::cout << "using custom arguments" << "\n";
        assert(argc == 3);
        char* dummy = NULL;
        memcpy(input_path, args[1], strlen(args[1]));
        memcpy(output_path, args[2], strlen(args[2]));
    }
    else {
        std::cout << "using default arguments" << "\n";
    }*/

    std::cout << "Arguments:" << "\n";
    std::cout << "\tin:\t" << input_path << "\n";
    std::cout << "\tout:\t" << output_path << "\n";
 
	std::cout << "level maker start" << std::endl;
    std::cout << "model at path " << LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_model.vox" << std::endl;

	const ogt_vox_scene* scene = load_vox_scene(LUSTRINE_EXPERIMENTS_DIRECTORY"/resources/level1_model.vox");
    assert(scene != nullptr);

    std::cout << "scene properties:" << std::endl;
    std::cout << "\tnum models:\t" << scene->num_models << std::endl;
    std::cout << "\tnum groups:\t" << scene->num_groups << std::endl;
    std::cout << "\tnum instances:\t" << scene->num_instances << std::endl;

    assert(scene->num_models > 0);

    const ogt_vox_model* model = scene->models[0];
    ogt_voxel_meshify_context ctx;
    memset(&ctx, 0, sizeof(ctx));
    ogt_mesh* mesh = ogt_mesh_from_paletted_voxels_simple(&ctx, model->voxel_data, model->size_x, model->size_y, model->size_z, (const ogt_mesh_rgba*)scene->palette.color);
    
    //assert(write_mesh_to_fbx(LUSTRINE_EXPERIMENTS_DIRECTORY"/voxel_to_obj/test.fbx", mesh, "mymesh"));

    std::unordered_map<glm::vec4, size_t> color_to_palette_index;
    glm::vec4* palette = (glm::vec4*) &scene->palette.color[0];

    for (uint32_t i = 0; i < 256; i++) {

        //color_to_palette_index[]

    }


    ogt_vox_destroy_scene(scene);
    std::cout << "level maker exit" << std::endl;
}