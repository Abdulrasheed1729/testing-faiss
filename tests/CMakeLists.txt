
add_executable(
main_test
scanner_test.cpp
kmer_utils_test.cpp
)

target_link_libraries(
main_test
PRIVATE lib_testing_faiss
GTest::gtest_main
)

target_compile_definitions(main_test PRIVATE
    TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/data"
)

set_target_properties(main_test PROPERTIES FOLDER tests)

include(GoogleTest)
gtest_discover_tests(main_test WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME main_test COMMAND main_test)

