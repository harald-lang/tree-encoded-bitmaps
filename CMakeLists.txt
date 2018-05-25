cmake_minimum_required(VERSION 3.8)
project(navin)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

add_subdirectory(../fastbit fastbit)
include_directories(../fastbit/src)

add_subdirectory(../CRoaring CRoaring)
include_directories(../CRoaring/include)

add_subdirectory(../sdsl-lite sdsl)
include_directories(../sdsl-lite/include)

add_subdirectory(../dtl dtl)
include_directories(../)

set(SOURCE_FILES
        index.hpp
        roaring_bitmap.hpp
        tree_mask_po.hpp
        two_state_markov_process.hpp
        static_stack.hpp
        util.hpp
        wah.hpp
        main.cpp
        bidirectional_mapping.hpp
        column_container.hpp
        column_imprints.hpp
        container_iterator.hpp
        file_parser.hpp
        h_psma.hpp
        histograms.hpp
        tests.hpp
        z_curve_runtime.hpp
        zCurve_compiletime.hpp
        )
add_executable(navin ${SOURCE_FILES})
target_link_libraries(navin fastbit roaring)
target_link_libraries(navin pthread)
target_link_libraries(navin m)
target_link_libraries(navin tbb)