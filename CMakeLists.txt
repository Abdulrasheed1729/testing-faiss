cmake_minimum_required(VERSION 4.1)
project(testing_faiss)

set(CMAKE_CXX_STANDARD 20)
find_package(faiss REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(testing_faiss main.cpp)
target_link_libraries(testing_faiss PRIVATE faiss)

