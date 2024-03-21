cmake_minimum_required(VERSION 3.16)

project(bgmtest LANGUAGES Fortran)
set(CMAKE_Fortran_FLAGS "-O2 -fconvert=big-endian")
#set(CMAKE_Fortran_PREPROCESS ON)

if(DEFINED W3LIBDIR)
link_directories(${W3LIBDIR})
else()
link_directories(lib)
endif()
if(DEFINED W3LIB)
set(_w3lib "${W3LIB}")
else()
set(_w3lib "w3_4")
endif()
add_subdirectory(common)

add_executable(addprtb addprtb.f90)
target_link_libraries(addprtb common)

add_executable(calcte calcte.f90)
target_link_libraries(calcte common)

include_directories(${CMAKE_BINARY_DIR}/common)
