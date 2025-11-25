cmake_minimum_required(VERSION 3.10)
project(JPEG C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compiler flags
add_compile_options(-Wall -Wextra)

# Core library source files
set(CORE_SOURCES
    core/block_process.c
    core/image_process.c
    core/streams.c
    core/constants.c
)

# Create a static library from the core sources
add_library(jpeg_core STATIC ${CORE_SOURCES})
target_include_directories(jpeg_core PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(jpeg_core m)

# Main executable (example/test program)
add_executable(jpeg_test main.c)
target_link_libraries(jpeg_test jpeg_core)

# Installation rules
install(TARGETS jpeg_core ARCHIVE DESTINATION lib)
install(FILES core/image_compression.h DESTINATION include/jpeg)
install(TARGETS jpeg_test RUNTIME DESTINATION bin)
