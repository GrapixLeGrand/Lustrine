project(particles)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma")


add_executable(wrapper wrapper.cpp)

target_link_libraries(wrapper LustrineWrapper)
