cmake_minimum_required(VERSION 3.12)
project(frame_utils)

# Set C++ standard (e.g., C++11)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the required packages
find_package(Eigen3 REQUIRED)

# Specify the locations of the header files
include_directories(include)

# # Collect source files
file(GLOB SOURCES "src/*.cpp")

# Create an executable target
#add_executable(convert_frames src/main.cpp)
add_executable(convert_frames ${SOURCES})

target_link_libraries(convert_frames Eigen3::Eigen)