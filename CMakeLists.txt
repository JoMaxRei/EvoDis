cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED on)

project(EvoDis)

find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

file(GLOB SOURCES "src/*.*")
add_executable(EvoDis ${SOURCES})

target_link_libraries(EvoDis GSL::gsl GSL::gslcblas)