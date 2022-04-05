project(wrapper)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -mfma")

#https://stackoverflow.com/questions/17511496/how-to-create-a-shared-library-with-cmake


add_executable(wrapper wrapper.cpp)

set(TARGET_INCLUDE_DIRS
	PUBLIC
    "../../src"
)

#target_include_directories(wrapper
#	${TARGET_INCLUDE_DIRS}
#)
#message(${CMAKE_BINARY_DIR}/LustrineWrapper.dll)
#add_library(LustrineWrapperImp SHARED IMPORTED GLOBAL)
# You can define two import-locations: one for debug and one for release.
#set_target_properties(LustrineWrapperImp PROPERTIES IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/LustrineWrapper.dll)
target_include_directories(
    wrapper
    PUBLIC "${Lustrine_SOURCE_DIR}/thirdparty/glm-0.9.9.8"
)
target_link_libraries(wrapper LustrineWrapper)
