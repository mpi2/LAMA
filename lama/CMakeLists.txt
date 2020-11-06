cmake_minimum_required(VERSION 3.11.0)

set(CMAKE_VERBOSE_MAKEFILE ON)

prjoect(itksnap_mesher)

find_package(itksnap)

if ( NOT itksnap_FOUND )
    message(STATUS "This project requires the ITKSNAP library, and will not be compiled")
    return()
endif()

# include helper file
include( ${ITKSNAP_USE_FILE} )

set(SOURCE_DIR "src/mesher")

include_directories(${SOURCE_DIR})

add_subdirectory(lib/pybind11)

pybind11_add_module(itksnap_mesher ${SOURCES} "${SOURCE_DIR}/itksnap_mesher.cpp")
