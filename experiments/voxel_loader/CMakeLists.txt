project(voxel_loader)

add_executable(voxel_loader voxel_loader.cpp)

target_link_libraries(voxel_loader Lustrine)

# Copy library DLL to target directory
add_custom_command(TARGET voxel_loader POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:voxel_loader>
        )