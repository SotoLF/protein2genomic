cmake_minimum_required(VERSION 3.16)
project(protein2genomic VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Find required packages
find_package(Threads REQUIRED)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/gtf_parser.cpp
    src/domain_mapper.cpp
    src/utils.cpp
)

# Create executable
add_executable(protein2genomic src/main.cpp ${SOURCES})

# Link libraries
target_link_libraries(protein2genomic Threads::Threads)

# Optional: Enable OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(protein2genomic OpenMP::OpenMP_CXX)
    target_compile_definitions(protein2genomic PRIVATE USE_OPENMP)
endif()

# Install
install(TARGETS protein2genomic DESTINATION bin)
