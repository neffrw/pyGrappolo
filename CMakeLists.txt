# Nitin A. Gawande, PNNL
# Oct 19, 2018

PROJECT(GRAPPOLO C CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.8...3.18)
add_subdirectory(pybind11)
set( TARGET grappolo)
pybind11_add_module(${TARGET} src/pyGrappolo.cpp)

# Include cmake modules
list( APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" )

# Project options
option( BUILD_SHARED_LIBS     
        "Build shared libraries instead of static libraries."
        OFF )
option( OPTION_SELF_CONTAINED 
        "Create a self-contained install with all dependencies."
        OFF )
option( OPTION_BUILD_TESTS
        "Build tests."
        ON )
option( OPTION_BUILD_DOCS  
        "Build documentation."
        OFF )
option( OPTION_BUILD_EXAMPLES 
        "Build examples."
        OFF )

set( CMAKE_CXX_STANDARD 11 )

#set(CMAKE_C_COMPILER "gcc-8")
#set(CMAKE_CXX_COMPILER "g++-8")
#set(CMAKE_C_COMPILER "/share/apps/gcc/7.1/bin/gcc")
#set(CMAKE_CXX_COMPILER "/share/apps/gcc/7.1/bin/g++")
# unset(CMAKE_SOURCE_DIR)
# set(CMAKE_SOURCE_DIR /grappolo)

set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -fopenmp -O3 -fPIC" )
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3 -fPIC" )
# set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -fopenmp -O0" )

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set varaibles using nomenclature used in Makefile
set( IOFOLDER grappolo/InputsOutput )
set( COFOLDER grappolo/BasicCommunitiesDetection )
set( UTFOLDER grappolo/Utility )
set( CLFOLDER grappolo/Coloring )
set( FSFOLDER grappolo/FullSyncOptimization )
set( DEFINSTR grappolo/DefineStructure )

# set(CMAKE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/grappolo)

set( IO_HEADERS
  ${CMAKE_SOURCE_DIR}/grappolo/DefineStructure
)

add_subdirectory( ${CMAKE_SOURCE_DIR}/grappolo/InputsOutput )
add_subdirectory( ${CMAKE_SOURCE_DIR}/grappolo/BasicCommunitiesDetection )
add_subdirectory( ${CMAKE_SOURCE_DIR}/grappolo/Utility )
add_subdirectory( ${CMAKE_SOURCE_DIR}/grappolo/Coloring )
add_subdirectory( ${CMAKE_SOURCE_DIR}/grappolo/FullSyncOptimization )

target_include_directories( ${TARGET} 
    PUBLIC ${IO_HEADERS}
)

target_link_libraries( ${TARGET} 
    PRIVATE
    inout  
    basic_cd
    util
    full_syn_opt
    coloring
)

target_compile_definitions(${TARGET}
    PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
