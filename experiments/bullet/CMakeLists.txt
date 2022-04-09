project(bullet)

add_executable(bullet bullet.cpp ../Utils.cpp)

if (WIN32)
    message("Lustrine::fluid: LevekGL path ${LevekGL_SOURCE_DIR}")
    set(GLEW_INCLUDE_DIRS "${LevekGL_SOURCE_DIR}/thirdparty/windows/glew-2.1.0/include")
    set(GLFW3_INCLUDE_DIRS "${LevekGL_SOURCE_DIR}/thirdparty/windows/glfw-3.3.6/include")

    set(GLEW_LIBRARIES "${LevekGL_SOURCE_DIR}/thirdparty/windows/glew-2.1.0/lib/Release/x64/glew32s.lib")
    set(GLFW3_LIBRARIES "${LevekGL_SOURCE_DIR}/thirdparty/windows/glfw-3.3.6/lib-vc2022/glfw3.lib")
    find_package(OpenGL REQUIRED)
    include_directories(${GLFW3_INCLUDE_DIRS})
    include_directories(${OPENGL_INCLUDE_DIR} ${GLEW_INCLUDE_DIRS})

    # Copy library DLL to target directory
    add_custom_command(TARGET bullet POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_FILE:Lustrine> $<TARGET_FILE_DIR:bullet>
            )
    add_custom_command(TARGET bullet POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${LevekGL_SOURCE_DIR}/thirdparty/windows/assimp/MD/assimp.dll $<TARGET_FILE_DIR:bullet>
            )
    add_custom_command(TARGET bullet POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            ${LevekGL_SOURCE_DIR}/thirdparty/windows/zlib-1.2.11/zlib.dll $<TARGET_FILE_DIR:bullet>
            )
endif (WIN32)

target_link_libraries(bullet Lustrine)