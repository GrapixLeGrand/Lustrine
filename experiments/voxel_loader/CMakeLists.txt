project(voxel_loader)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma")

add_executable(voxel_loader voxel_loader.cpp)

target_link_libraries(voxel_loader Lustrine)
