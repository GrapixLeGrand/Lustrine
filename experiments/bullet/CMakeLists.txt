project(bullet)

add_executable(bullet bullet.cpp ../Utils.cpp)

target_link_libraries(bullet Lustrine)

# Copy library DLL to target directory
#add_custom_command(TARGET voxel_loader POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy_if_different
#        $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:bullet>
#        )