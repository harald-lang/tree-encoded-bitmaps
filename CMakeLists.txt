cmake_minimum_required(VERSION 3.5)
project(navin)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native -Wno-sign-compare")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native -mavx512bw -mavx512f")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=core-avx2")

# HACK: This gives access to the private members of the boost::dynamic_bitset.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS")

# Enable ASAN in debug builds.
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

include_directories(src)
include_directories(repro/src)

add_subdirectory(../fastbit Fastbit) # Fastbit: which includes WAH
include_directories(SYSTEM ../fastbit/src)

# Roaring Bitmaps (deprecated since we switched to the unity build)
#add_subdirectory(../CRoaring CRoaring)
#include_directories(../CRoaring/include)

#include_directories(../sdsl-lite/include)

add_subdirectory(../dtl dtl)
include_directories(SYSTEM ../dtl/src)
#include_directories(dtl_storage/src)
include_directories(../) #TODO remove

include_directories(./)

add_subdirectory(thirdparty/gtest)
include_directories(thirdparty/gtest/include lib/googletest)
include_directories(thirdparty)


##===----------------------------------------------------------------------===##
# Expose Git commit ID.
#
# Inspired from
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/
# and
# https://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automat
# ically-print-out-its-git-version-hash
##===----------------------------------------------------------------------===##

# Get the current working branch.
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch.
execute_process(
        COMMAND git describe --abbrev=7 --dirty=+ --always --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the current date and time.
execute_process(
        COMMAND date +%Y-%m-%d
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BUILD_DATE
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process(
        COMMAND date +%H:%M:%S
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE BUILD_TIME
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

#add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
#add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")
#add_definitions("-DBUILD_DATE=${BUILD_DATE}")
#add_definitions("-DBUILD_TIME=${BUILD_TIME}")

configure_file(
        ${CMAKE_SOURCE_DIR}/src/version.h.in
        ${CMAKE_BINARY_DIR}/generated/version.h
)
include_directories(${CMAKE_BINARY_DIR}/generated)
##===----------------------------------------------------------------------===##

set(SQLITE_SOURCE_FILES
        thirdparty/sqlite/sqlite3.c
        )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSQLITE_OMIT_LOAD_EXTENSION")

set(OTHER_SOURCE_FILES
        thirdparty/perfevent/PerfEvent.hpp
        thirdparty/types/Types.cpp
        thirdparty/types/Types.hpp
        thirdparty/roaring/roaring.c
        )

set(BENCHMARK_SOURCE_FILES
        ${SQLITE_SOURCE_FILES}
        ${OTHER_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/performance/common.hpp
        experiments/performance/common_bitwise.hpp
        experiments/util/bitmap_db.hpp
        experiments/util/bitmap_db.cpp
        experiments/util/bitmap_types.hpp
        experiments/util/config.hpp
        experiments/util/gen.hpp
        experiments/util/gen.cpp
        experiments/util/params.hpp
        experiments/util/prep_data.cpp
        experiments/util/prep_data.hpp
        experiments/util/seq_db.hpp
        experiments/util/seq_db.cpp
        experiments/util/threading.hpp
        experiments/util/types.hpp
        )

##===----------------------------------------------------------------------===##
# Test cases.
##===----------------------------------------------------------------------===##
set(TEST_FILES
        ${BENCHMARK_SOURCE_FILES}
        test/dtl/bitmap/util/bit_buffer_avx2_test.cpp
        test/dtl/bitmap/util/bit_buffer_avx512_test.cpp
        test/dtl/bitmap/util/bit_buffer_test.cpp
        test/dtl/bitmap/util/bitmap_fun_test.cpp
        test/dtl/bitmap/util/rank_test.cpp
        test/dtl/bitmap/api_types.hpp
        test/dtl/bitmap/api_encode_decode_test.cpp
        test/dtl/bitmap/api_random_access_test.cpp
        test/dtl/bitmap/api_run_iterator_test.cpp
        test/dtl/bitmap/api_run_iterator_skip_test.cpp
        test/dtl/bitmap/api_bitwise_operation_test.cpp
        test/dtl/bitmap/bitwise_operations_helper.hpp
        test/dtl/bitmap/diff_test.cpp
        test/dtl/bitmap/part_diff_test.cpp
        test/dtl/bitmap/update_test.cpp
        test/dtl/bitmap/teb_scan_util_test.cpp
        )
add_executable(tester ${TEST_FILES})
target_link_libraries(tester gtest gtest_main fastbit dl)

##===----------------------------------------------------------------------===##
# Sources
##===----------------------------------------------------------------------===##
set(UTIL_SOURCE_FILES
        repro/src/dtl/bitmap/util/base64.hpp
        repro/src/dtl/bitmap/util/markov_process.hpp
        repro/src/dtl/bitmap/util/random.hpp
        repro/src/dtl/bitmap/util/two_state_markov_process.hpp
        src/dtl/bitmap/util/convert.hpp
        )

set(TEB_SOURCE_FILES
        src/dtl/bitmap.hpp
        src/dtl/bitmap/diff/diff.hpp
        src/dtl/bitmap/diff/merge.hpp
        src/dtl/bitmap/diff/merge_teb.hpp
        src/dtl/bitmap/part/part.hpp
        src/dtl/bitmap/part/part_run.hpp
        src/dtl/bitmap/part/part_updirect.hpp
        src/dtl/bitmap/part/part_upforward.hpp
        src/dtl/bitmap/util/binary_tree_structure.hpp
        src/dtl/bitmap/util/bit_buffer.hpp
        src/dtl/bitmap/util/bit_buffer_avx2.hpp
        src/dtl/bitmap/util/bit_buffer_avx512.hpp
        src/dtl/bitmap/util/bitmap_fun.hpp
        src/dtl/bitmap/util/bitmap_tree.hpp
        src/dtl/bitmap/util/bitmap_view.hpp
        src/dtl/bitmap/util/mutable_bitmap_tree.hpp
        src/dtl/bitmap/util/plain_bitmap.hpp
        src/dtl/bitmap/util/rank1_naive.hpp
        src/dtl/bitmap/util/rank1.hpp
        src/dtl/bitmap/util/rank1_logic_surf.hpp
        src/dtl/bitmap/bitwise_operations.hpp
        src/dtl/bitmap/iterator.hpp
        src/dtl/bitmap/teb.hpp
        src/dtl/bitmap/teb_builder.hpp
        src/dtl/bitmap/teb_flat.hpp
        src/dtl/bitmap/teb_iter.hpp
        src/dtl/bitmap/teb_wrapper.hpp
        src/dtl/bitmap/teb_scan_iter.hpp
        src/dtl/bitmap/teb_scan_util.hpp
        src/dtl/bitmap/teb_types.hpp
        src/dtl/iterator.hpp
        src/dtl/static_stack.hpp
        src/dtl/static_stack2.hpp
        )

# Our competitors
set(COMPETITORS_SOURCE_FILES
        # Boost dynamic bitset
        repro/src/dtl/bitmap/dynamic_bitmap.hpp

        # Roaring
        repro/src/dtl/bitmap/dynamic_roaring_bitmap.hpp

        # WAH (FastBit)
        repro/src/dtl/bitmap/dynamic_wah.hpp

        # Position list
        repro/src/dtl/bitmap/position_list.hpp
        repro/src/dtl/bitmap/partitioned_position_list.hpp

        # Range list
        repro/src/dtl/bitmap/range_list.hpp
        repro/src/dtl/bitmap/partitioned_range_list.hpp
        )

set(SOURCE_FILES
        ${UTIL_SOURCE_FILES}
        ${TEB_SOURCE_FILES}
        ${COMPETITORS_SOURCE_FILES}
        )

##===----------------------------------------------------------------------===##
# Tools, benchmarks, and experiments.
##===----------------------------------------------------------------------===##
# Space optimizations. - Compares the different optimizations levels.
set(EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/main_compression_optimizations.cpp
        )
add_executable(ex_compression_opt ${EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES})
target_link_libraries(ex_compression_opt fastbit pthread dl)

# TODO add description
set(EXPERIMENT_COMPRESSION_OPT_METRICS_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/main_compression_optimizations_metrics.cpp
        )
add_executable(ex_compression_opt_metrics ${EXPERIMENT_COMPRESSION_OPT_METRICS_SOURCE_FILES})
target_link_libraries(ex_compression_opt_metrics fastbit pthread dl)

# Compression, skyline (markov)
set(EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline.cpp
        )
add_executable(ex_compression_skyline ${EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline fastbit pthread dl)

# Compression, skyline (markov, fine grained)
set(EXPERIMENT_COMPRESSION_SKYLINE_FINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline_fine.cpp
        )
add_executable(ex_compression_skyline_fine ${EXPERIMENT_COMPRESSION_SKYLINE_FINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline_fine fastbit pthread dl)

# Compression, varying d (uniform)
set(EXPERIMENT_COMPRESSION_UNIFORM_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_uniform.cpp
        )
add_executable(ex_compression_uniform ${EXPERIMENT_COMPRESSION_UNIFORM_SOURCE_FILES})
target_link_libraries(ex_compression_uniform fastbit pthread dl)

# Compression, varying f (markov)
set(EXPERIMENT_COMPRESSION_VARYING_F_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_varying_f.cpp
        )
add_executable(ex_compression_varying_f ${EXPERIMENT_COMPRESSION_VARYING_F_SOURCE_FILES})
target_link_libraries(ex_compression_varying_f fastbit pthread dl)

# Compression, varying d (markov)
set(EXPERIMENT_COMPRESSION_VARYING_D_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_varying_d.cpp
        )
add_executable(ex_compression_varying_d ${EXPERIMENT_COMPRESSION_VARYING_D_SOURCE_FILES})
target_link_libraries(ex_compression_varying_d fastbit pthread dl)

# Compression, real data sets.
set(REAL_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/real/real.cpp
        )
add_executable(real ${REAL_SOURCE_FILES})
target_link_libraries(real fastbit pthread dl boost_system boost_filesystem jemalloc)

# copy the data set to the build directory
file(COPY experiments/real/real-roaring-dataset DESTINATION ${PROJECT_BINARY_DIR})

# Performance, varying f (markov)
set(EXPERIMENT_PERFORMANCE_VARYING_F_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_varying_f.cpp
        )
add_executable(ex_performance_varying_f ${EXPERIMENT_PERFORMANCE_VARYING_F_SOURCE_FILES})
target_link_libraries(ex_performance_varying_f fastbit pthread dl)

# Performance, skyline
set(EXPERIMENT_PERFORMANCE_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_performance_skyline.cpp
        )
add_executable(ex_performance_skyline ${EXPERIMENT_PERFORMANCE_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_performance_skyline fastbit pthread dl)

# Performance, intersect
set(EXPERIMENT_PERFORMANCE_INTERSECT_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_intersect.cpp
        )
add_executable(ex_performance_intersect ${EXPERIMENT_PERFORMANCE_INTERSECT_SOURCE_FILES})
target_link_libraries(ex_performance_intersect fastbit pthread dl)

set(EXPERIMENT_MICROBENCHMARK_DOWNWARDS_NAV_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_microbenchmark_downwards_navigation.cpp
        )
add_executable(ex_microbenchmark_downwards_nav ${EXPERIMENT_MICROBENCHMARK_DOWNWARDS_NAV_SOURCE_FILES})
target_link_libraries(ex_microbenchmark_downwards_nav fastbit pthread dl)

set(EXPERIMENT_MICROBENCHMARK_UPWARDS_NAV_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_microbenchmark_upwards_navigation.cpp
        )
add_executable(ex_microbenchmark_upwards_nav ${EXPERIMENT_MICROBENCHMARK_UPWARDS_NAV_SOURCE_FILES})
target_link_libraries(ex_microbenchmark_upwards_nav fastbit pthread dl)

# Index compression
set(EXPERIMENT_INDEX_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/index/main_index_compression.cpp
        )
add_executable(ex_index_compression ${EXPERIMENT_INDEX_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_index_compression fastbit pthread dl)

##===----------------------------------------------------------------------===##
include(CMakeListsLocal.cmake OPTIONAL)
