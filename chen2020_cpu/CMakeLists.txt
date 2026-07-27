cmake_minimum_required(VERSION 3.16)
project(pbcd_rgod LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(CHEN2020_ENABLE_OPENMP "Build the PBCD implementation with OpenMP" ON)

add_library(chen2020_core STATIC
    src/TNM.cpp
    src/TNM_Net.cpp
    src/TNM_Algorithm.cpp
    src/TNM_EnumIO.cpp
    src/TNM_utility.cpp
    src/PBCD_Algorithm.cpp
)

target_compile_definitions(chen2020_core PUBLIC TNM_SHUTMYSQL)
target_include_directories(chen2020_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/randlib
)

if(CHEN2020_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(chen2020_core PUBLIC OpenMP::OpenMP_CXX)
endif()

if(UNIX)
    target_link_libraries(chen2020_core PUBLIC m)
endif()

add_executable(chen2020_ue src/driver.cpp)
target_link_libraries(chen2020_ue PRIVATE chen2020_core)

include(CTest)
if(BUILD_TESTING)
    add_executable(pbcd_block_test tests/pbcd_block_test.cpp)
    target_link_libraries(pbcd_block_test PRIVATE chen2020_core)
    add_test(NAME pbcd_block_test COMMAND pbcd_block_test)

    add_executable(pbcd_rgod_test tests/pbcd_rgod_test.cpp)
    target_link_libraries(pbcd_rgod_test PRIVATE chen2020_core)
    add_test(NAME pbcd_rgod_test COMMAND pbcd_rgod_test)

    add_executable(tntp_sparse_origin_test tests/tntp_sparse_origin_test.cpp)
    target_link_libraries(tntp_sparse_origin_test PRIVATE chen2020_core)
    add_test(NAME tntp_sparse_origin_test COMMAND tntp_sparse_origin_test)
endif()

set_target_properties(chen2020_ue PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
