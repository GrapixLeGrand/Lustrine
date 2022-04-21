project(performance)

add_executable(performance main.cpp "tsc_x86.h")

if (WIN32)
    # Copy library DLL to target directory
    add_custom_command(TARGET performance POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:bullet>
            )
    add_custom_command(TARGET performance POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${LevekGL_SOURCE_DIR}/thirdparty/windows/assimp/MD/assimp.dll $<TARGET_FILE_DIR:bullet>
            )
    add_custom_command(TARGET performance POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${LevekGL_SOURCE_DIR}/thirdparty/windows/zlib-1.2.11/zlib.dll $<TARGET_FILE_DIR:bullet>
            )
endif (WIN32)

target_link_libraries(performance Lustrine)