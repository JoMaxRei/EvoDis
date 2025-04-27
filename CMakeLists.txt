cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)

project(EvoDis)

find_package(GSL REQUIRED)

include_directories(${GSL_INCLUDE_DIRS})
link_directories(${GSL_LIBRARY_DIRS})

file(GLOB SOURCES "src/*.cpp*" "src/*.cc")
add_executable(EvoDis ${SOURCES})


target_link_libraries(EvoDis GSL::gsl GSL::gslcblas)
get_target_property(GSL_DLL_PATH GSL::gsl IMPORTED_LOCATION)
add_custom_command(TARGET EvoDis POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${GSL_DLL_PATH} $<TARGET_FILE_DIR:EvoDis>)