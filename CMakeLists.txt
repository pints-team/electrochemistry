cmake_minimum_required(VERSION 3.1)

project(electrochemistry)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1 )

set(include_dirs "")
set(lib_dir "lib")

# Eigen
list(APPEND include_dirs ${lib_dir}/eigen3)

# Pybind11
add_subdirectory(${lib_dir}/pybind11)

# Boost 
list(APPEND include_dirs ${lib_dir}/boost)

set(source_dir src/electrochemistry)

set(source_files 
    ${source_dir}/e_implicit_exponential_mesh.cpp
    ${source_dir}/seq_electron_transfer3_explicit.cpp
    )
set(header_files 
    ${source_dir}/e_implicit_exponential_mesh.hpp
    ${source_dir}/seq_electron_transfer3_explicit.hpp
    ${source_dir}/utilities.hpp
    )

include_directories(${source_dir} ${include_dirs})

pybind11_add_module(electrochemistry ${source_files} ${header_files} "${source_dir}/wrapper.cpp")

