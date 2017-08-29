project(TriGA_CPP)
cmake_minimum_required(VERSION 3.9.0)

# set the output folder where the project executable will be
# set(CMAKE_BINARY_DIR ${CMAKESOURCE_DIR}/BIN)
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
# set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories(include)

# add subfolders
# ------------------------
file(GLOB SOURCES "lib/*.cpp")
# add_subdirectory (gnuplot)
# add_subdirectory (IO_files)

add_executable(TriGA_CPP ${SOURCES}