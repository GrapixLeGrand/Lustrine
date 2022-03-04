project(particles)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma")

if (UNIX)
	#
	#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
	#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
endif (UNIX)

add_executable(fluid fluid.cpp ../Utils.cpp)

target_link_libraries(fluid Lustrine)
