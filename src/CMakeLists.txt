# file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${testing_faiss_SOURCE_DIR}/include/*.hpp")

add_library(lib_testing_faiss
        faiss_utils.cpp
        # io.cpp
        scanner.cpp
)

target_include_directories(lib_testing_faiss PUBLIC
        ${testing_faiss_SOURCE_DIR}/include
)

target_link_libraries(lib_testing_faiss
        faiss
        OpenMP::OpenMP_CXX
)
