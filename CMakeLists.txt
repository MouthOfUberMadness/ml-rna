cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set the project name
project(ml-rna VERSION 0.1)

# include Bio++ libraries
find_package (bpp-core3 1.0.0 REQUIRED)
find_package (bpp-seq3 1.0.0 REQUIRED)
find_package (bpp-phyl3 1.0.0 REQUIRED)

# add the executable
add_executable(ml-rna src/main.cpp src/analysis.cpp)

# link to libraries
target_link_libraries(ml-rna ${BPP_LIBS_SHARED})