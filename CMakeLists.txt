cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

project(UKBioCoin LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(UKC_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

find_package(Boost REQUIRED COMPONENTS program_options)
find_package(OpenMP REQUIRED COMPONENTS CXX)

set(EIGEN3_INCLUDE_DIR "" CACHE PATH "Directory containing Eigen/Dense")
if(NOT EIGEN3_INCLUDE_DIR)
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Dense
        HINTS ENV EIGEN3_INCLUDE_DIR /home/jingcheng/fsl/include/eigen3
        PATH_SUFFIXES eigen3)
endif()
if(NOT EIGEN3_INCLUDE_DIR)
    message(FATAL_ERROR "Eigen was not found. Set -DEIGEN3_INCLUDE_DIR=/path/to/eigen3")
endif()

add_executable(UKBioCoin main.cpp)
target_include_directories(UKBioCoin PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${EIGEN3_INCLUDE_DIR}")
target_link_libraries(UKBioCoin PRIVATE
    Boost::program_options
    OpenMP::OpenMP_CXX)
target_compile_options(UKBioCoin PRIVATE -Wall -Wextra -Wpedantic)

if(UKC_ENABLE_SANITIZERS)
    target_compile_options(UKBioCoin PRIVATE
        -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
    target_link_options(UKBioCoin PRIVATE -fsanitize=address,undefined)
endif()

include(CTest)
if(BUILD_TESTING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(ukbiocoin_cdf_smoke tests/test_cdf.cpp)
    target_include_directories(ukbiocoin_cdf_smoke PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ukbiocoin_cdf_smoke COMMAND ukbiocoin_cdf_smoke)

    add_test(
        NAME ukbiocoin_regression
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ukbiocoin.py"
                "$<TARGET_FILE:UKBioCoin>")
    set_tests_properties(ukbiocoin_regression PROPERTIES
        ENVIRONMENT "OMP_THREAD_LIMIT=16"
        TIMEOUT 120)
endif()
