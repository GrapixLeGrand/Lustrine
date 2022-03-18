project(particles)

add_executable(fluid fluid.cpp ../Utils.cpp)

target_link_libraries(fluid Lustrine)
