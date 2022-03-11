project(particles)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma")

if (UNIX)
	#
	#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
	#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif (UNIX)

if (WIN32)
set(GLEW_INCLUDE_DIRS "../../thirdparty/LevekGL/thirdparty/windows/glew-2.1.0/include")
set(GLFW3_INCLUDE_DIRS "../../thridparty/LevekGL/thirdparty/windows/glfw-3.3.6/include")
#set(GLEW_LIBRARIES "${LevekGL_SOURCE_DIR}/thirdparty/windows/glew-2.1.0/lib/Release/x64/glew32s.lib")
#set(GLFW3_LIBRARIES "${LevekGL_SOURCE_DIR}/thirdparty/windows/glfw-3.3.6/lib-vc2022/glfw3.lib")
include_directories(GLEW_INCLUDE_DIRS)
include_directories(GLFW3_INCLUDE_DIRS)
#target_include_directories(fluid ${GLEW_INCLUDE_DIRS} ${GLFW3_INCLUDE_DIR
endif(WIN32)

add_executable(fluid fluid.cpp ../Utils.cpp)

target_link_libraries(fluid Lustrine)
