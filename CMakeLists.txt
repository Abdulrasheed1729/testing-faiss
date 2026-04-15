cmake_minimum_required(VERSION 4.1)
project(testing_faiss)

set(CMAKE_CXX_STANDARD 20)
find_package(faiss REQUIRED)
find_package(OpenMP REQUIRED)
find_package(GTest)

add_executable(testing_faiss main.cpp)
target_link_libraries(testing_faiss PRIVATE faiss)

enable_testing()

add_executable(
parser_test
parser_test.cpp
)

target_link_libraries(
parser_test
GTest::gtest_main
)

target_compile_definitions(parser_test PRIVATE
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

include(GoogleTest)
gtest_discover_tests(parser_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
