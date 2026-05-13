cmake_minimum_required(VERSION 3.16)
project(prot2exon VERSION 1.0.0 LANGUAGES CXX)

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
    src/output_writer.cpp
    src/utils.cpp
)

# Create executable
add_executable(prot2exon src/main.cpp ${SOURCES})

# Link libraries
target_link_libraries(prot2exon Threads::Threads)

# Optional: Enable OpenMP for parallel processing
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(prot2exon OpenMP::OpenMP_CXX)
    target_compile_definitions(prot2exon PRIVATE USE_OPENMP)
endif()

# Install
install(TARGETS prot2exon DESTINATION bin)
