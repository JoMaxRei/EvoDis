cmake_minimum_required(VERSION 3.11)

project(EvoDis)

find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

file(GLOB SOURCES "src/*.cpp")
add_executable(EvoDis ${SOURCES})

target_link_libraries(EvoDis GSL::gsl GSL::gslcblas)