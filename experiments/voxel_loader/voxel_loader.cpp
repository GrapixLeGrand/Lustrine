#include <iostream>
#include <fstream>
#include <vector>

#include "VoxelLoader.hpp"

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: voxel_loader VOX_FILE\n";
        return 1;
    }

    std::ifstream fh(argv[1]);
    if (fh) {
        fh.seekg(std::ios_base::end);
        int64_t size = fh.tellg();
        fh.seekg(std::ios_base::beg);
        std::vector<uint8_t> data(size);
        fh.read((char*) data.data(), size);
        std::cout << Lustrine::read_vox_scene_json(data.data(), size);
        return 0;
    } else {
        std::cerr << "Failed to open the file!\n";
        return 2;
    }
}