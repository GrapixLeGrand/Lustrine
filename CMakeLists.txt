cmake_minimum_required(VERSION 3.13)
project(Lustrine VERSION 0.1.0)
cmake_policy(SET CMP0079 NEW)
cmake_policy(SET CMP0111 NEW)  # error on missing imported targets
include(CTest)
include(GenerateExportHeader)
enable_testing()

if (UNIX)
    add_compile_definitions(PLATFORM_UNIX)
endif()

if (WIN32)
    add_compile_definitions(PLATFORM_WINDOWS)
endif()


#option(ENV{LEVEKGL_USE_INTERNAL_GLM} "Include paths to external GLM library" ON)
#option(ENV{LEVEKGL_USE_INTERNAL_GLM} OFF)
#option(ENV{LEVEKGL_EXTERNAL_GLM} "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8/glm/glm.hpp")

set(ENV{LEVEKGL_USE_INTERNAL_GLM_REQUESTED})
set(ENV{LEVEKGL_EXTERNAL_GLM} "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8")

add_subdirectory(thirdparty)

link_directories(${Lustrine_SOURCE_DIR}/src)

set(LUSTRINE_SOURCES
    "src/Lustrine.cpp"
    "src/JsonWriter.cpp"
    "src/Kernels.cpp"
    "src/Neighbors.cpp"
    "src/VoxelLoader.cpp"
    "src/BulletPhysics.cpp"
)

set(LUSTRINEWRAPPER_SOURCES
    "src/Lustrine.cpp"
    "src/JsonWriter.cpp"
    "src/Kernels.cpp"
    "src/Neighbors.cpp"
    "src/LustrineWrapper.cpp"
    "src/BulletPhysics.cpp"
    "src/VoxelLoader.cpp"
)


add_library(Lustrine SHARED ${LUSTRINE_SOURCES})
generate_export_header(Lustrine
        BASE_NAME Lustrine
        EXPORT_MACRO_NAME LUSTRINE_EXPORT
        EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/Lustrine_Export.h
        STATIC_DEFINE LUSTRINE_BUILT_AS_STATIC
)
target_include_directories(Lustrine PUBLIC ${CMAKE_BINARY_DIR}/exports)

add_library(LustrineWrapper SHARED ${LUSTRINEWRAPPER_SOURCES})
#add_library(LustrineWrapper STATIC ${LUSTRINEWRAPPER_SOURCES})
generate_export_header(LustrineWrapper
        BASE_NAME LustrineWrapper
        EXPORT_MACRO_NAME LUSTRINE_WRAPPER_EXPORT
        EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/LustrineWrapper_Export.h
        STATIC_DEFINE LUSTRINE_WRAPPER_BUILT_AS_STATIC
)
set_target_properties(
        LustrineWrapper PROPERTIES
        PUBLIC_HEADER src/LustrineWrapper.hpp
)
target_include_directories(LustrineWrapper PUBLIC ${CMAKE_BINARY_DIR}/exports)


target_link_libraries(Lustrine LevekGL)


if (WIN32)
    set(BULLET_INCLUDE_DIR
            "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/bullet3/src"
            )
    set(Bullet3LinearMathLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/LinearMath.lib")
    set(Bullet3CommonLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/Bullet3Common.lib")
    set(BulletDynamicsLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletDynamics.lib")
    set(BulletInverseDynamicsLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletInverseDynamics.lib")
    set(BulletCollisionLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletCollision.lib")
    set(BulletSoftBodyLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletSoftBody.lib")

    add_library(Bullet3LinearMath STATIC IMPORTED)
    set_property(TARGET Bullet3LinearMath PROPERTY IMPORTED_LOCATION ${Bullet3LinearMathLibPath})
    list (APPEND BULLET_LIBRARIES LinearMath)

    add_library(Bullet3Common STATIC IMPORTED)
    set_property(TARGET Bullet3Common PROPERTY IMPORTED_LOCATION ${Bullet3CommonLibPath})
    list (APPEND BULLET_LIBRARIES Bullet3Common)

    add_library(Bullet3InverseDynamics STATIC IMPORTED)
    set_property(TARGET Bullet3InverseDynamics PROPERTY IMPORTED_LOCATION ${BulletInverseDynamicsLibPath})
    list (APPEND BULLET_LIBRARIES BulletInverseDynamics)

    add_library(BulletCollision STATIC IMPORTED)
    set_property(TARGET BulletCollision PROPERTY IMPORTED_LOCATION ${BulletCollisionLibPath})
    list (APPEND BULLET_LIBRARIES BulletCollision)

    add_library(BulletDynamics STATIC IMPORTED)
    set_property(TARGET BulletDynamics PROPERTY IMPORTED_LOCATION ${BulletDynamicsLibPath})
    list (APPEND BULLET_LIBRARIES BulletDynamics)

    add_library(BulletSoftBody STATIC IMPORTED)
    set_property(TARGET BulletSoftBody PROPERTY IMPORTED_LOCATION ${BulletSoftBodyLibPath})
    list (APPEND BULLET_LIBRARIES BulletSoftBody)

    message("linear math ${Bullet3LinearMathLibPath}")
    message("common ${Bullet3CommonLibPath}")
    message("Bullet3InverseDynamics ${BulletInverseDynamicsLibPath}")
    message("BulletCollision ${BulletCollisionLibPath}")
    message("BulletDynamics ${BulletDynamicsLibPath}")
    message("Bullet3SoftBody ${BulletSoftBodyLibPath}")
else()
    find_package(Bullet REQUIRED)

    message("bullet libraries: ${BULLET_LIBRARIES}")
    message("bullet include dir: ${BULLET_INCLUDE_DIR}")
endif()

target_link_libraries(Lustrine ${BULLET_LIBRARIES})
target_link_libraries(LustrineWrapper ${BULLET_LIBRARIES})
target_include_directories(Lustrine PUBLIC ${BULLET_INCLUDE_DIR})
target_include_directories(LustrineWrapper PUBLIC ${BULLET_INCLUDE_DIR})

set(LUSTRINE_TARGET_INCLUDE_DIRS
    PUBLIC "src"
    PRIVATE "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8"#"thirdparty/glm-0.9.9.8"
    PUBLIC "thirdparty/LevekGL"   
)

set(LUSTRINEWRAPPER_TARGET_INCLUDE_DIRS
    PUBLIC "src"
    PRIVATE "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8"#"thirdparty/glm-0.9.9.8"   
)

target_include_directories(Lustrine ${LUSTRINE_TARGET_INCLUDE_DIRS})
target_include_directories(LustrineWrapper ${LUSTRINEWRAPPER_TARGET_INCLUDE_DIRS})

add_compile_definitions(LUSTRINE_EXPERIMENTS_DIRECTORY="${Lustrine_SOURCE_DIR}/experiments")

add_subdirectory(experiments)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
