cmake_minimum_required(VERSION 3.10)
project(hspist3 LANGUAGES C CXX)

# Add the GLFW include directory
include_directories("/opt/homebrew/Cellar/glfw/3.4/include")

# Add the GLFW library directory
link_directories("/opt/homebrew/Cellar/glfw/3.4/lib")

# Add the source files
file(GLOB SOURCES "src/*.c")

# Create the executable
add_executable(hspist3 ${SOURCES})

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Find packages
find_package(OpenGL REQUIRED)
find_package(SDL2 REQUIRED)

# Link libraries
target_link_libraries(hspist3 SDL2::SDL2 OpenGL::GL glfw)