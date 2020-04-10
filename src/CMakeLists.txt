cmake_minimum_required (VERSION 3.0)
set(CMAKE_CXX_STANDARD 14)
project(strukern)
add_library(my_python3.7m SHARED IMPORTED)
set_target_properties(my_python3.7m PROPERTIES
        IMPORTED_LOCATION "/home/mrjoeybux/programs/miniconda3/envs/mypy/lib/libpython3.7m.so"
        INTERFACE_INCLUDE_DIRECTORIES "/home/mrjoeybux/programs/miniconda3/envs/mypy/include/python3.7m"
        )
find_package(PythonLibs REQUIRED)
find_package(pybind11 REQUIRED)
set(EIGEN3_INCLUDE_DIR /home/mrjoeybux/programs/eigen/include/eigen3/)
include_directories(${EIGEN3_INCLUDE_DIR})
pybind11_add_module(strukern binding/strukern.cpp)