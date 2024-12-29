cmake_minimum_required(VERSION 3.11)

project(EvoDis)

file(GLOB SOURCES "src/*.cpp")
add_executable(EvoDis ${SOURCES})
