project(performance)

add_executable(performance main.cpp)

if (WIN32)
    # Copy library DLL to target directory
    add_custom_command(TARGET performance POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:performance>
            )
endif (WIN32)

target_link_libraries(performance Lustrine)