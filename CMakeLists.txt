cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

project(UKBioCoin VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(UKC_ENABLE_SANITIZERS "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)

find_package(Boost REQUIRED COMPONENTS charconv filesystem json program_options)
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

add_library(ukbiocoin_nss STATIC nss_v2.cpp)
target_include_directories(ukbiocoin_nss PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(ukbiocoin_nss PUBLIC Boost::filesystem Boost::json)
target_compile_options(ukbiocoin_nss PRIVATE -Wall -Wextra -Wpedantic)

add_executable(UKBioCoin main.cpp)
target_include_directories(UKBioCoin PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${EIGEN3_INCLUDE_DIR}")
target_link_libraries(UKBioCoin PRIVATE
    Boost::charconv
    Boost::program_options
    ukbiocoin_nss
    OpenMP::OpenMP_CXX)
target_compile_options(UKBioCoin PRIVATE -Wall -Wextra -Wpedantic)

add_executable(UKBioCoinNSSAdapter nss_adapter.cpp)
target_include_directories(UKBioCoinNSSAdapter PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(UKBioCoinNSSAdapter PRIVATE
    Boost::program_options ukbiocoin_nss OpenMP::OpenMP_CXX)
target_compile_options(UKBioCoinNSSAdapter PRIVATE -Wall -Wextra -Wpedantic)

if(UKC_ENABLE_SANITIZERS)
    foreach(target UKBioCoin UKBioCoinNSSAdapter ukbiocoin_nss)
        target_compile_options(${target} PRIVATE
            -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined)
        target_link_options(${target} PRIVATE -fsanitize=address,undefined)
    endforeach()
endif()

include(CTest)
if(BUILD_TESTING)
    if(DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/bin/python")
        set(Python3_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python" CACHE FILEPATH
            "Python interpreter used by tests")
    endif()
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
    add_executable(ukbiocoin_cdf_smoke tests/test_cdf.cpp)
    target_include_directories(ukbiocoin_cdf_smoke PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ukbiocoin_cdf_smoke COMMAND ukbiocoin_cdf_smoke)

    add_test(
        NAME ukbiocoin_regression
        COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_ukbiocoin.py"
                "$<TARGET_FILE:UKBioCoin>"
                "$<TARGET_FILE:UKBioCoinNSSAdapter>")
    set_tests_properties(ukbiocoin_regression PROPERTIES
        ENVIRONMENT "OMP_THREAD_LIMIT=16"
        TIMEOUT 120)
endif()
