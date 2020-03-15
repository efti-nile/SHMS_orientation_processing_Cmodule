cmake_minimum_required(VERSION 3.14)
project(SHMS_orientation_processing_Cmodule C)

set(CMAKE_C_STANDARD 99)

# These files aren't gonna be changed so I use GLOB
FILE(GLOB LINEAR_ALGEBRA_SOURCES ${CMAKE_SOURCE_DIR}/LinearAlgebra/*.c)
FILE(GLOB LAPACK_SOURCES
        ${CMAKE_SOURCE_DIR}/Lapack/Blas/*.c
        ${CMAKE_SOURCE_DIR}/Lapack/Scr/*.c
        ${CMAKE_SOURCE_DIR}/Lapack/F2c/*.c
        ${CMAKE_SOURCE_DIR}/Lapack/Install/*.c
        )
FILE(GLOB QPOASES_SOURCES ${CMAKE_SOURCE_DIR}/qpOASES/Source/*.c)

add_executable(SHMS_orientation_processing_Cmodule main.c orientation_processing.c orientation_processing.h gauss_newton_algo.c gauss_newton_algo.h)
target_link_libraries(SHMS_orientation_processing_Cmodule m)
add_library(linear_algebra STATIC ${LINEAR_ALGEBRA_SOURCES} ${LAPACK_SOURCES} ${QPOASES_SOURCES})
target_link_libraries(SHMS_orientation_processing_Cmodule linear_algebra)
