cmake_minimum_required(VERSION 3.8)
project(navin)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=corei7")

include_directories(src)


add_subdirectory(../fastbit Fastbit) # Fastbit includes the WAH bitmap compression
include_directories(../fastbit/src)

add_subdirectory(../CRoaring CRoaring)
include_directories(../CRoaring/include)

#include_directories(../sdsl-lite/include)

add_subdirectory(../dtl dtl)
include_directories(../dtl/src)
include_directories(dtl_storage/src)
include_directories(../) #TODO remove

include_directories(./)

add_subdirectory(thirdparty/gtest)
include_directories(thirdparty/gtest/include lib/googletest)
include_directories(thirdparty)


set(SQLITE_SOURCE_FILES
        thirdparty/sqlite/sqlite3.c
        )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSQLITE_OMIT_LOAD_EXTENSION")

set(BENCHMARK_SOURCE_FILES
        ${SQLITE_SOURCE_FILES}
        experiments/util/database.hpp
        experiments/util/database.cpp
        experiments/util/gen.hpp
        )

##===----------------------------------------------------------------------===##
# Test cases.
##===----------------------------------------------------------------------===##
set(TEST_FILES
        test/dtl/bitmap/api_encode_decode_test.cpp
        test/dtl/bitmap/api_encode_decode_using_iterator_test.cpp
        test/dtl/bitmap/teb_test.cpp
        test/dtl/bitmap/api_range_iterator_skip_test.cpp
#        test/dtl/bitmap/api_set_operation_test.cpp
#        test/dtl/bitmap/static_api_test.cpp
#        bitmap_view_test.cpp
#        test/dtl/bitmap/iterate_test.cpp
        )
add_executable(tester ${TEST_FILES})
target_link_libraries(tester gtest gtest_main fastbit roaring)

# TODO removed
#set(TEST_ITER_FILES
#        test/dtl/bitmap/iter_tester.cpp)
#add_executable(tester_iter ${TEST_ITER_FILES})
#target_link_libraries(tester_iter gtest gtest_main)
##===----------------------------------------------------------------------===##

##===----------------------------------------------------------------------===##
# Sources
##===----------------------------------------------------------------------===##
set(TEB_SOURCE_FILES
        src/dtl/bitmap.hpp
        src/dtl/bitmap/dynamic_partitioned_tree_mask.hpp
#        src/dtl/bitmap/dynamic_tree_mask_lo.hpp
        src/dtl/bitmap/util/base64.hpp
        src/dtl/bitmap/util/bitmap_tree.hpp
        src/dtl/bitmap/util/binary_tree_structure.hpp
        src/dtl/bitmap/util/convert.hpp
        src/dtl/bitmap/util/rank1.hpp
        src/dtl/bitmap/util/random.hpp
        src/dtl/bitmap/util/two_state_markov_process.hpp
        src/dtl/bitmap/teb.hpp
        src/dtl/bitmap/tree_mask_util.hpp
        src/dtl/static_stack.hpp
        )

# Our competitors
set(COMPETITORS_SOURCE_FILES
        src/dtl/bitmap/dynamic_bitmap.hpp
        src/dtl/bitmap/dynamic_roaring_bitmap.hpp
        src/dtl/bitmap/dynamic_wah.hpp
        src/dtl/bitmap/partitioned_range_list.hpp
        src/dtl/bitmap/partitioned_position_list.hpp
        src/dtl/bitmap/position_list.hpp
        src/dtl/bitmap/range_list.hpp
        )

set(SOURCE_FILES
        ${TEB_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        #bitmap_view.hpp
#        index.hpp
#        util.hpp
#        src/dtl/storage/bidirectional_mapping.hpp
#        src/dtl/storage/column_container.hpp
#        src/dtl/storage/column_imprints.hpp
#        src/dtl/storage/container_iterator.hpp
#        src/dtl/storage/file_parser.hpp
#        src/dtl/storage/h_psma.hpp
#        src/dtl/storage/histograms.hpp
#        test/dtl/storage/tests.hpp
#        src/dtl/storage/z_curve_runtime.hpp
#        src/dtl/storage/zCurve_compiletime.hpp
        )

#add_executable(navin main.cpp ${SOURCE_FILES})
#target_link_libraries(navin fastbit roaring)
#target_link_libraries(navin pthread)
#target_link_libraries(navin m)
#target_link_libraries(navin tbb)



##===----------------------------------------------------------------------===##
# Tools, benchmarks, and experiments.
##===----------------------------------------------------------------------===##
#set(GENERATOR_SOURCE_FILES
#        ${SQLITE_SOURCE_FILES}
#        experiments/util/database.hpp
#        experiments/util/database.cpp
#        tools/gen/gen_bitmaps.cpp
#        )
#add_executable(gen_bitmaps ${GENERATOR_SOURCE_FILES})
#target_link_libraries(gen_bitmaps pthread dl)

# TBD
set(EXPERIMENT_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression.cpp
        )
add_executable(ex_compression ${EXPERIMENT_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_compression fastbit roaring)

# Space optimizations. - Compares the different optimizations levels.
set(EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_optimizations.cpp
        )
add_executable(ex_compression_opt ${EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES})
target_link_libraries(ex_compression_opt fastbit roaring)

# Skyline
set(EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline.cpp
        )
add_executable(ex_compression_skyline ${EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline fastbit roaring pthread dl)


# copy shell scripts
configure_file(experiments/compression/compression.sh compression.sh COPYONLY)


set(EXPERIMENT_PERFORMANCE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_performance.cpp
        )
add_executable(ex_performance ${EXPERIMENT_PERFORMANCE_SOURCE_FILES})
target_link_libraries(ex_performance fastbit roaring)


set(EXPERIMENT_LOSSY_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        experiments/compression/main_lossycompression.cpp
        )
add_executable(ex_lossycompression ${EXPERIMENT_LOSSY_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_lossycompression fastbit roaring)
##===----------------------------------------------------------------------===##
