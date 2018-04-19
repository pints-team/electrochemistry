cmake_minimum_required(VERSION 2.8.3)
project(electrochemistry)

set(include_dirs "")
set(lib_dir "lib")

# Eigen
list(APPEND ${lib_dir}/eigen3)

# Pybind11
add_subdirectory(${lib_dir}/pybind11)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1 )
list(APPEND CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++14")

set(source_dir src/electrochemistry)

set(source_files 
    ${source_dir}/e_implicit_exponential_mesh.cpp
    ${source_dir}/src/electrochemistry/seq_electron_transfer3_explicit.cpp
    )
set(header_files 
    ${source_dir}/e_implicit_exponential_mesh.hpp
    ${source_dir}/seq_electron_transfer3_explicit.hpp
    ${source_dir}/utilities.hpp
    )

include_directories(${source_dir} ${lib_dir})

pybind11_add_module(electrochemistry ${source_files} ${header_files} "${source_dir}/wrapper.cpp")

