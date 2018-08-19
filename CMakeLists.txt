cmake_minimum_required(VERSION 3.8)
project(navin)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native")

add_subdirectory(../fastbit Fastbit) # Fastbit includes the WAH bitmap compression
include_directories(../fastbit/src)

add_subdirectory(../CRoaring CRoaring)
include_directories(../CRoaring/include)

include_directories(../sdsl-lite/include)

add_subdirectory(../dtl dtl)
include_directories(../dtl/src)
include_directories(dtl_storage/src)
include_directories(../) #TODO remove

add_subdirectory(thirdparty/gtest)
include_directories(thirdparty/gtest/include lib/googletest)

set(TEST_FILES
        api_test.cpp
#        bitmap_view_test.cpp
        iterate_1fills_test.cpp
        )
add_executable(tester ${TEST_FILES})
target_link_libraries(tester gtest gtest_main fastbit roaring)

# Treemask sources
set(TREE_MASK_SOURCE_FILES
        dynamic_bitmap.hpp
        dynamic_partitioned_tree_mask.hpp
        dynamic_tree.hpp
        dynamic_tree_mask_lo.hpp
        partitioned_tree_mask.hpp
        rank.hpp
        tree_mask_lo.hpp
        tree_mask_po.hpp
        tree_mask_util.hpp
        )

# Our third-party competitors
set(COMPETITORS_SOURCE_FILES
        bitmap.hpp # std::bitset
        roaring_bitmap.hpp
        wah.hpp
        )

set(SOURCE_FILES
        ${TREE_MASK_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        bitmap_view.hpp
        index.hpp
        two_state_markov_process.hpp
        static_stack.hpp
        util.hpp
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

add_executable(navin main.cpp ${SOURCE_FILES})
target_link_libraries(navin fastbit roaring)
target_link_libraries(navin pthread)
target_link_libraries(navin m)
target_link_libraries(navin tbb)


set(EXPERIMENT_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        ${TREE_MASK_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression.cpp
        )

add_executable(ex_compression ${EXPERIMENT_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_compression fastbit roaring)


set(EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${TREE_MASK_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline.cpp
        )
add_executable(ex_compression_skyline ${EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline fastbit roaring)


# copy shell scripts
configure_file(experiments/compression/compression.sh compression.sh COPYONLY)

set(EXPERIMENT_PERFORMANCE_SOURCE_FILES
        ${SOURCE_FILES}
        ${TREE_MASK_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        experiments/performance/main_performance.cpp
        )

add_executable(ex_performance ${EXPERIMENT_PERFORMANCE_SOURCE_FILES})
target_link_libraries(ex_performance fastbit roaring)
