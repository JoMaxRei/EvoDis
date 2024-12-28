cmake_minimum_required(VERSION 3.10)

project(EvoDis)

file(GLOB SOURCES "src/*.cpp")
add_executable(EvoDis ${SOURCES})