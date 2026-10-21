cmake_minimum_required(VERSION 3.24)
set(PROJECT_LANGUAGES "CXX")
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_FORTRAN_STANDARD 90)
set(PROJECT_FAMILY_NAME "Tennsix")
set(PROJECT_VERSION "3.1.2")
set(PROJECT_VERSION_NAME "MultiSquid")
set(PROJECT_VERSION_MONTH "February")
set(PROJECT_VERSION_YEAR "2023")
set(PROJECT_WORKING_RDIRS "_T" "_I")
# Release, Debug
set(PROJECT_BUILD_TYPE Debug)
set(ACCEL "") # "Cuda", "Metal", "" (none)
set(ACCEL_DEPENDENCIES "") # cuBLAS, cuRAND, ...
set(ACCEL_INSTALL_DIR "/usr/local") # /usr/local, $ENV{CUDA_HOME}, $ENV{METAL_HOME}, ...
set(DEPENDENCY_LIB_TYPE "s") # shared (s) or static (a)
set(DGCALC_DEPENDENCY_TEST "")
set(DGCALC_DEPENDENCY_SRC ${DGCALC_DEPENDENCY_TEST} "Base" "Kyte" "Types" "Scalar")
set(DGCALC_DEPENDENCY_LIB "")
set(GROUPS "THIRD")
set(THIRD_DEPENDENCY_TEST "Catch2")
set(THIRD_DEPENDENCY_SRC ${THIRD_DEPENDENCY_TEST})
set(THIRD_DEPENDENCY_LIB "")

# set(LOCAL_ROOT "")
if (LOCAL_ROOT)
    include(${LOCAL_ROOT}/src/dgcalc.1.cmake)
else()
    execute_process(COMMAND dgcalc source 1 OUTPUT_VARIABLE SOURCE1_CODE)
    cmake_language(EVAL CODE "${SOURCE1_CODE}")
endif()
set(PROJECT_NAME "${PROJECT_FAMILY_NAME}${PROJECT_VERSION}${PROJECT_VERSION_NAME}_${PROJECT_VERSION_MONTH}${PROJECT_VERSION_YEAR}")
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES ${PROJECT_LANGUAGES})
if (LOCAL_ROOT)
    include(${LOCAL_ROOT}/src/dgcalc.2.cmake)
else()
    execute_process(COMMAND dgcalc source 2 OUTPUT_VARIABLE SOURCE2_CODE)
    cmake_language(EVAL CODE "${SOURCE2_CODE}")
endif()






