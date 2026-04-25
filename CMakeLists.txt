cmake_minimum_required(VERSION 4.1)
project(testing_faiss)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

find_package(faiss REQUIRED)
find_package(OpenMP REQUIRED)


add_executable(testing_faiss main.cpp)
add_subdirectory(src)
target_link_libraries(testing_faiss 
	PRIVATE faiss 
	PRIVATE lib_testing_faiss
)


# Testing
find_package(GTest)
enable_testing()
add_subdirectory(tests)

# add_executable(
# scanner_test
# scanner_test.cpp
# )
#
# target_link_libraries(
# scanner_test
# GTest::gtest_main
# )
#
# target_compile_definitions(scanner_test PRIVATE
#     TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
# )
#
# include(GoogleTest)
# gtest_discover_tests(scanner_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
#
# add_executable(
# kmer_utils_test
# kmer_utils_test.cpp
# )
#
# target_link_libraries(
# kmer_utils_test
# GTest::gtest_main
# )
#
# gtest_discover_tests(kmer_utils_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
