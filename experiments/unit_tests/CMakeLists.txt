project(unit_testing)

add_executable(unit_testing main.cpp)

if (WIN32)
    # Copy library DLL to target directory
    add_custom_command(TARGET unit_testing POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:unit_testing>
            )
endif (WIN32)

target_link_libraries(unit_testing Lustrine)