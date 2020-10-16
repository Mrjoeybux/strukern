cmake_minimum_required (VERSION 3.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
#set(SOURCE_DIR "./")
#include_directories({SOURCE_DIR})
message("${CMAKE_SOURCE_DIR}")
project(strukern)
find_package(PythonLibs REQUIRED)
add_subdirectory(lib/pybind11)
find_package(ZLIB)
#find_package(JPEG REQUIRED)
find_package (Eigen3 REQUIRED NO_MODULE)

#set(GEDLIB /home/Mrjoeybux/programs/gedlib/)
set(GEDLIB ${CMAKE_SOURCE_DIR}/lib/gedlib/)
set(GEDLIB_BOOST ${GEDLIB}/ext/boost.1.69.0/)
set(GEDLIB_EIGEN ${GEDLIB}/ext/eigen.3.3.4/Eigen/)
set(GEDLIB_NOMAD ${GEDLIB}/ext/nomad.3.8.1/src/)
set(GEDLIB_NOMAD_SGTE ${GEDLIB}/ext/nomad.3.8.1/ext/sgtelib/src/)
set(GEDLIB_LSAPE ${GEDLIB}/ext/lsape.5/include/)
set(GEDLIB_LIBSVM ${GEDLIB}/ext/libsvm.3.22/)
set(GEDLIB_FANN ${GEDLIB}/ext/fann.2.2.0/include/)
#find_package(dlib REQUIRED)
set(DLIB ${CMAKE_SOURCE_DIR}/lib/dlib/)
#set(DLIB /home/Mrjoeybux/coding/libs/dlib-19.16)

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
set(strukern_src main.cpp)
pybind11_add_module(strukern binding/strukern.cpp)

set(DLIB_LINK ${DLIB}/build/dlib/libdlib.so)
set(NOMAD_LINK ${GEDLIB}/ext/nomad.3.8.1/lib/libnomad.so)
set(LIBSVM_LINK ${GEDLIB}/ext/libsvm.3.22/libsvm.so)
set(FANN_LINK ${GEDLIB}/ext/fann.2.2.0/lib/libdoublefann.so)
set(GEDLIB_GXL_LINK ${GEDLIB}/lib/libgxlgedlib.so)
target_link_libraries(strukern PRIVATE ${GEDLIB_GXL_LINK} ${ZLIB_LIBRARIES} ${NOMAD_LINK} ${LIBSVM_LINK} ${FANN_LINK} ${JPEG_LIBRARIES} ${DLIB_LINK})# domain kernels)

