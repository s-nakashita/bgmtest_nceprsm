cmake_minimum_required(VERSION 3.16)

file(GLOB_RECURSE sources *.f90)
add_library(common STATIC ${sources})
target_link_libraries(common ${_w3lib})
#target_compile_options(common PRIVATE -fno-range-check)
