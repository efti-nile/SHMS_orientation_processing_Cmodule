cmake_minimum_required(VERSION 3.14)
project(SHMS_orientation_processing_Cmodule C)

set(CMAKE_C_STANDARD 99)

add_executable(SHMS_orientation_processing_Cmodule main.c orientation_processing.c orientation_processing.h gauss_newton_algo.c gauss_newton_algo.h linear_algebra.c linear_algebra.h)
target_link_libraries(SHMS_orientation_processing_Cmodule m)
