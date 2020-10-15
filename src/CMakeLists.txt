cmake_minimum_required (VERSION 3.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS
#        -DCMAKE_BUILD_TYPE=Release
#        -DPYTHON_EXECUTABLE=$HOME/miniconda/envs/krein/bin/python3.6
#        -DPYTHON_LIBRARY=$HOME/miniconda/envs/krein/lib/libpython3.6m.so
#        -DPYTHON_INCLUDE_DIR=$HOME/miniconda/envs/krein/include/python3.6m)
project(strukern)
#add_library(my_python3.7m SHARED IMPORTED)
#set(PYTHON_INCLUDE_DIRS /home/Mrjoeybux/programs/miniconda3/include/python3.6m/)
#find_package(PythonLibs REQUIRED)
#include(FindPythonAnaconda)
set(PYTHONLIBS_FOUND TRUE)
set(PYTHON_MODULE_EXTENSION .so)
set(PYTHON_MODULE_PREFIX "")
#find_package(pybind11 REQUIRED)
add_subdirectory(pybind11)
find_package(ZLIB)
find_package(JPEG REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)
#find_package(PythonLibs REQUIRED)
set(GEDLIB /home/Mrjoeybux/programs/gedlib/)
set(GEDLIB_BOOST ${GEDLIB}/ext/boost.1.69.0/)
set(GEDLIB_EIGEN ${GEDLIB}/ext/eigen.3.3.4/Eigen/)
set(GEDLIB_NOMAD ${GEDLIB}/ext/nomad.3.8.1/src/)
set(GEDLIB_NOMAD_SGTE ${GEDLIB}/ext/nomad.3.8.1/ext/sgtelib/src/)
set(GEDLIB_LSAPE ${GEDLIB}/ext/lsape.5/include/)
set(GEDLIB_LIBSVM ${GEDLIB}/ext/libsvm.3.22/)
set(GEDLIB_FANN ${GEDLIB}/ext/fann.2.2.0/include/)
set(MY_EIGEN /home/mrjoeybux/programs/eigen/include/eigen3/)
set(DLIB /home/Mrjoeybux/coding/libs/dlib-19.16)
#set(PYTHON_INCLUDE_PATH "-DPYTHON_EXECUTABLE")
set(PYTHON_INCLUDE_PATH /home/Mrjoeybux/programs/miniconda3/include/python3.6m/)
set(PYTHON_LIBRARIES "/home/Mrjoeybux/programs/miniconda3/lib/libpython3.6m.so")
set(PYTHON_INCLUDE_DIR /home/Mrjoeybux/programs/miniconda3/include/python3.6m/)
include_directories(${GEDLIB} 
                ${GEDLIB_BOOST} 
                ${GEDLIB_EIGEN} 
                ${GEDLIB_NOMAD} 
                ${GEDLIB_NOMAD_SGTE} 
                ${GEDLIB_LSAPE} 
                ${GEDLIB_LIBSVM} 
                ${GEDLIB_FANN}
                ${PYTHON_INCLUDE_DIR}
                ${JPEG_INCLUDE_DIR}
                ${EIGEN3_INCLUDE_DIR}
                ${DLIB}
                )
                #${MY_EIGEN})
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "INCLUDED DIRECTORIES:")
foreach(dir ${dirs})
        message(STATUS "${dir}")
endforeach()


set(NOMAD_LIB ${GEDLIB}/ext/nomad.3.8.1/lib/)
set(LIBSVM_LIB ${GEDLIB}/ext/libsvm.3.22/)
set(FANN_LIB ${GEDLIB}/ext/fann.2.2.0/lib/)
set(GEDLIB_GXL_LIB ${GEDLIB}/lib/)

link_directories(${GEDLIB_GXL_LIB} ${NOMAD_LIB} ${LIBSVM_LIB} ${FANN_LIB} ${PYTHON_LIBRARIES})
#add_subdirectory(domain)
#add_subdirectory(kernels)
set(strukern_src main.cpp)
#add_executable(strukern ${strukern_src})

#set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "ged_test")
pybind11_add_module(strukern binding/strukern.cpp)

set(DLIB_LINK ${DLIB}/build/dlib/libdlib.so)
set(NOMAD_LINK ${GEDLIB}/ext/nomad.3.8.1/lib/libnomad.so)
set(LIBSVM_LINK ${GEDLIB}/ext/libsvm.3.22/libsvm.so)
set(FANN_LINK ${GEDLIB}/ext/fann.2.2.0/lib/libdoublefann.so)
set(GEDLIB_GXL_LINK ${GEDLIB}/lib/libgxlgedlib.so)
message(${GEDLIB_GXL_LINK})
#target_link_libraries(strukern PRIVATE ${ZLIB_LIBRARIES} ${JPEG_LIBRARIES} ${GEDLIB_GXL_LINK})
target_link_libraries(strukern PRIVATE ${GEDLIB_GXL_LINK} ${ZLIB_LIBRARIES} ${NOMAD_LINK} ${LIBSVM_LINK} ${FANN_LINK} ${JPEG_LIBRARIES} ${DLIB_LINK})# domain kernels)

