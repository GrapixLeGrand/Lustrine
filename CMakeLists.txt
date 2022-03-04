cmake_minimum_required(VERSION 3.13)
project(Lustrine VERSION 0.1.0)
cmake_policy(SET CMP0079 NEW)
include(CTest)
enable_testing()

#option(ENV{LEVEKGL_USE_INTERNAL_GLM} "Include paths to external GLM library" ON)
#option(ENV{LEVEKGL_USE_INTERNAL_GLM} OFF)
#option(ENV{LEVEKGL_EXTERNAL_GLM} "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8/glm/glm.hpp")

set(ENV{LEVEKGL_USE_INTERNAL_GLM_REQUESTED})
set(ENV{LEVEKGL_EXTERNAL_GLM} "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8/glm/glm.hpp")

add_subdirectory(thirdparty)

link_directories(${Lustrine_SOURCE_DIR}/src)

set(LUSTRINE_SOURCES
    "src/Lustrine.cpp"
    "src/Kernels.cpp"
    "src/Neighbors.cpp"
)

add_library(Lustrine STATIC ${LUSTRINE_SOURCES})
target_link_libraries(Lustrine LevekGL)

set(LUSTRINE_TARGET_INCLUDE_DIRS
    PUBLIC "src"
    PRIVATE "thirdparty/glm-0.9.9.8"
    PUBLIC "thirdparty/LevekGL"   
)

target_include_directories(Lustrine ${LUSTRINE_TARGET_INCLUDE_DIRS})

add_compile_definitions(LUSTRINE_EXPERIMENTS_DIRECTORY="${Lustrine_SOURCE_DIR}/experiments")

add_subdirectory(experiments)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
