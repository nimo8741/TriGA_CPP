project(TriGA_CPP)
cmake_minimum_required(VERSION 3.9.0)

# set the output folder where the project executable will be
set(ROOT_DIR ${PROJECT_SOURCE_DIR})
set(EXECUTABLE_OUTPUT_PATH ${ROOT_DIR}/build)
set(SOURCE_DIR ${PROJECT_SOURCE_DIR}/lib)
MESSAGE ("Source path: ${ROOT_DIR}")

include_directories(Eigen_lib)
include_directories(include)

# add subfolders
# ------------------------
add_subdirectory (lib)

add_executable(mesh_gen ${SOURCE_DIR}/mesh_main.cpp ${SOURCE_DIR}/Mesh_Part1_funcs.cpp ${SOURCE_DIR}/Mesh_Part2_funcs.cpp ${SOURCE_DIR}/Mesh_Part3_funcs.cpp ${SOURCE_DIR}/Mesh_Part4_funcs.cpp)