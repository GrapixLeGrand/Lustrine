cmake_minimum_required(VERSION 3.13)
project(Lustrine VERSION 0.1.0)
cmake_policy(SET CMP0079 NEW)
include(CTest)
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
    "src/Kernels.cpp"
    "src/Neighbors.cpp"
    "src/VoxelLoader.cpp"
    "src/BulletPhysics.cpp"
)

set(LUSTRINEWRAPPER_SOURCES
    "src/Lustrine.cpp"
    "src/Kernels.cpp"
    "src/Neighbors.cpp"
    "src/LustrineWrapper.cpp"
    "src/BulletPhysics.cpp"
)


add_library(Lustrine SHARED ${LUSTRINE_SOURCES})
add_library(LustrineWrapper SHARED ${LUSTRINEWRAPPER_SOURCES})
#add_library(LustrineWrapper STATIC ${LUSTRINEWRAPPER_SOURCES})
set_target_properties(
    LustrineWrapper PROPERTIES
    PUBLIC_HEADER src/LustrineWrapper.hpp
)

target_link_libraries(Lustrine LevekGL)


set(Bullet3IncludeDirs 
    PUBLIC "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/bullet3/src"
)
if (WIN32)
#thirdparty\Bullet3\windows
 
#    BulletDynamics
#    BulletInverseDynamics
#    BulletCollision
#    BulletSoftBody

set(Bullet3LinearMathLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/LinearMath.lib")
set(Bullet3CommonLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/Bullet3Common.lib")
set(BulletDynamicsLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletDynamics.lib")
set(BulletInverseDynamicsLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletInverseDynamics.lib")
set(BulletCollisionLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletCollision.lib")
set(BulletSoftBodyLibPath "${Lustrine_SOURCE_DIR}/thirdparty/Bullet3/windows/MD/BulletSoftBody.lib")

endif()

if (UNIX)
#TODO set the static lib for bullet here for unix
endif()

target_include_directories(
    Lustrine
    ${Bullet3IncludeDirs}
)

target_include_directories(
    LustrineWrapper
    ${Bullet3IncludeDirs}
)
#include_directories(${Bullet3IncludeDirs})




add_library(Bullet3LinearMath STATIC IMPORTED)
set_property(TARGET Bullet3LinearMath PROPERTY IMPORTED_LOCATION ${Bullet3LinearMathLibPath})
target_link_libraries(Lustrine Bullet3LinearMath)
target_link_libraries(LustrineWrapper Bullet3LinearMath)

add_library(Bullet3Common STATIC IMPORTED)
set_property(TARGET Bullet3Common PROPERTY IMPORTED_LOCATION ${Bullet3CommonLibPath})
target_link_libraries(Lustrine Bullet3Common)
target_link_libraries(LustrineWrapper Bullet3Common)

add_library(Bullet3InverseDynamics STATIC IMPORTED)
set_property(TARGET Bullet3InverseDynamics PROPERTY IMPORTED_LOCATION ${BulletInverseDynamicsLibPath})
target_link_libraries(Lustrine Bullet3InverseDynamics)
target_link_libraries(LustrineWrapper Bullet3InverseDynamics)

add_library(BulletCollision STATIC IMPORTED)
set_property(TARGET BulletCollision PROPERTY IMPORTED_LOCATION ${BulletCollisionLibPath})
target_link_libraries(Lustrine BulletCollision)
target_link_libraries(LustrineWrapper BulletCollision)

add_library(BulletDynamics STATIC IMPORTED)
set_property(TARGET BulletDynamics PROPERTY IMPORTED_LOCATION ${BulletDynamicsLibPath})
target_link_libraries(Lustrine BulletDynamics)
target_link_libraries(LustrineWrapper BulletDynamics)


add_library(Bullet3SoftBody STATIC IMPORTED)
set_property(TARGET Bullet3SoftBody PROPERTY IMPORTED_LOCATION ${BulletSoftBodyLibPath})
target_link_libraries(Lustrine Bullet3SoftBody)
target_link_libraries(LustrineWrapper Bullet3SoftBody)


message("linear math ${Bullet3LinearMathLibPath}")
message("common ${Bullet3CommonLibPath}")
message("Bullet3InverseDynamics ${BulletInverseDynamicsLibPath}")
message("BulletCollision ${BulletCollisionLibPath}")
message("BulletDynamics ${BulletDynamicsLibPath}")
message("Bullet3SoftBody ${BulletSoftBodyLibPath}")

#list (APPEND BULLET_DEFINITIONS ${BULLET_DOUBLE_DEF})
#list (APPEND BULLET_LIBRARIES LinearMath)
#list (APPEND BULLET_LIBRARIES Bullet3Common)
#list (APPEND BULLET_LIBRARIES BulletInverseDynamics)
#list (APPEND BULLET_LIBRARIES BulletCollision)
#list (APPEND BULLET_LIBRARIES BulletDynamics)
#list (APPEND BULLET_LIBRARIES BulletSoftBody)

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
