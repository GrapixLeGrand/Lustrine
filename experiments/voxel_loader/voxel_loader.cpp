#include <iostream>
#include <fstream>
#include <vector>

#include "VoxelLoader.hpp"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: voxel_loader VOX_FILE\n";
        return 1;
    }

    std::ifstream fh(argv[1], std::ios::binary | std::ios::ate);
    if (fh) {
        int64_t size = fh.tellg();
        fh.seekg(std::ios_base::beg);
        std::cerr << "Reading file (" << size << " bytes)" << std::endl;
        std::vector<uint8_t> data(size);
        fh.read((char*) data.data(), size);
        std::string json = Lustrine::read_vox_scene_json(data.data(), size);
        std::cerr << "Writing to console..." << std::endl;
        std::cout << json;
        return 0;
    } else {
        std::cerr << "Failed to open the file!\n";
        return 2;
    }
}