cmake_minimum_required(VERSION 3.12)
project(Detcresolution)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 

# Find ROOT
find_package(ROOT REQUIRED)

# Include ROOT headers
include_directories(${ROOT_INCLUDE_DIRS})

# Add source files
set(SOURCES DigiValidation.cc DigiChecker.cc)

# Link ROOT libraries
add_executable(DigiValidation ${SOURCES})
target_link_libraries(DigiValidation ${ROOT_LIBRARIES})
