cmake_minimum_required(VERSION 3.10)

# Project name and language
project(particle_assembler LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the executable (replace "main.cpp" with your actual filename)
add_executable(${PROJECT_NAME} source/particle_assembler.cc)
