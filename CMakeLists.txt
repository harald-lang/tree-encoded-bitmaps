cmake_minimum_required(VERSION 3.5)
project(navin)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native -Wno-sign-compare")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=native -mavx512bw -mavx512f")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -march=core-avx2")

# Hack: This gives us access to the private members of the dynamic_bitset.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS")

# Enable ASAN in debug builds.
#set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
#set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

include_directories(src)


add_subdirectory(../fastbit Fastbit) # Fastbit includes the WAH bitmap compression
include_directories(../fastbit/src)

# Roaring Bitmaps
#add_subdirectory(../CRoaring CRoaring)
#include_directories(../CRoaring/include)

#include_directories(../sdsl-lite/include)

add_subdirectory(../dtl dtl)
include_directories(../dtl/src)
include_directories(dtl_storage/src)
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
        experiments/util/bitmap_types.hpp
        experiments/util/bitmap_db.hpp
        experiments/util/bitmap_db.cpp
        experiments/util/config.hpp
        experiments/util/gen.hpp
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
#        test/dtl/bitmap/api_encode_decode_test.cpp /* DEPRECATED */
        test/dtl/bitmap/api_encode_decode_using_iterator_test.cpp
        test/dtl/bitmap/api_random_access_test.cpp
        test/dtl/bitmap/bit_buffer_test.cpp
        test/dtl/bitmap/bit_buffer_avx2_test.cpp
        test/dtl/bitmap/bit_buffer_avx512_test.cpp
        test/dtl/bitmap/bitmap_fun_test.cpp
        test/dtl/bitmap/teb_test.cpp
        test/dtl/bitmap/teb_scan_util_test.cpp
        test/dtl/bitmap/api_range_iterator_skip_test.cpp
        test/dtl/bitmap/api_set_operation_test.cpp
        test/dtl/bitmap/util/bitmap_tree_test.cpp
        test/dtl/bitmap/util/bitwise_operations.hpp
#        test/dtl/bitmap/static_api_test.cpp
#        test/dtl/bitmap/iterate_test.cpp
        )
add_executable(tester ${TEST_FILES})
target_link_libraries(tester gtest gtest_main fastbit dl)

set(VALIDATOR_FILES
        ${BENCHMARK_SOURCE_FILES}
        test/dtl/bitmap/encode_decode_pregenerated_bitmaps_test.cpp
        )
add_executable(validator ${VALIDATOR_FILES})
target_link_libraries(validator gtest gtest_main fastbit dl)

# To be removed
set(TEST_FILES
        ${BENCHMARK_SOURCE_FILES}
        test/dtl/bitmap/sandbox_test.cpp
        )
add_executable(sandbox_test ${TEST_FILES})
target_link_libraries(sandbox_test gtest gtest_main fastbit dl)

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
        src/dtl/bitmap/util/bit_buffer.hpp
        src/dtl/bitmap/util/bit_buffer_avx2.hpp
        src/dtl/bitmap/util/bit_buffer_avx512.hpp
        src/dtl/bitmap/util/plain_bitmap.hpp
        src/dtl/bitmap/util/bitmap_fun.hpp
        src/dtl/bitmap/util/bitmap_tree.hpp
        src/dtl/bitmap/util/bitmap_view.hpp
        src/dtl/bitmap/util/binary_tree_structure.hpp
        src/dtl/bitmap/util/convert.hpp
        src/dtl/bitmap/util/rank1.hpp
        src/dtl/bitmap/util/rank1_interceptor.hpp
        src/dtl/bitmap/util/rank1_naive.hpp
        src/dtl/bitmap/util/rank1_surf.hpp
        src/dtl/bitmap/util/rank1_surf_cached.hpp
        src/dtl/bitmap/util/rank1_super_fast.hpp
        src/dtl/bitmap/util/random.hpp
        src/dtl/bitmap/util/markov_process.hpp
        src/dtl/bitmap/util/two_state_markov_process.hpp
        src/dtl/bitmap/bitwise_operations.hpp
        src/dtl/bitmap/teb.hpp
        src/dtl/bitmap/teb_flat.hpp
        src/dtl/bitmap/teb_wrapper.hpp
        src/dtl/bitmap/teb_scan.hpp
        src/dtl/bitmap/teb_scan_util.hpp
        src/dtl/bitmap/tree_mask_util.hpp
        src/dtl/iterator.hpp
        src/dtl/static_stack.hpp
        src/dtl/static_stack2.hpp
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

# Space optimizations. - Compares the different optimizations levels.
set(EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/main_compression_optimizations.cpp
        )
add_executable(ex_compression_opt ${EXPERIMENT_COMPRESSION_OPT_SOURCE_FILES})
target_link_libraries(ex_compression_opt fastbit pthread dl)

set(EXPERIMENT_COMPRESSION_OPT_MATRICS_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/main_compression_optimizations_metrics.cpp
        )
add_executable(ex_compression_opt_metrics ${EXPERIMENT_COMPRESSION_OPT_MATRICS_SOURCE_FILES})
target_link_libraries(ex_compression_opt_metrics fastbit pthread dl)

# Skyline (markov)
set(EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline.cpp
        )
add_executable(ex_compression_skyline ${EXPERIMENT_COMPRESSION_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline fastbit pthread dl)

# Skyline (markov, fine grained)
set(EXPERIMENT_COMPRESSION_SKYLINE_FINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/compression/main_compression_skyline_fine.cpp
        )
add_executable(ex_compression_skyline_fine ${EXPERIMENT_COMPRESSION_SKYLINE_FINE_SOURCE_FILES})
target_link_libraries(ex_compression_skyline_fine fastbit pthread dl)

# Compression (uniform)
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

# Performance, varying f (markov)
set(EXPERIMENT_PERFORMANCE_VARYING_F_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_varying_f.cpp
        )
add_executable(ex_performance_varying_f ${EXPERIMENT_PERFORMANCE_VARYING_F_SOURCE_FILES})
target_link_libraries(ex_performance_varying_f fastbit pthread dl)

# Performance, construction
set(EXPERIMENT_PERFORMANCE_CONSTRUCTION_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_construction.cpp
        )
add_executable(ex_performance_construction ${EXPERIMENT_PERFORMANCE_CONSTRUCTION_SOURCE_FILES})
target_link_libraries(ex_performance_construction fastbit pthread dl)

# Performance, intersect
set(EXPERIMENT_PERFORMANCE_INTERSECT_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_intersect.cpp
        )
add_executable(ex_performance_intersect ${EXPERIMENT_PERFORMANCE_INTERSECT_SOURCE_FILES})
target_link_libraries(ex_performance_intersect fastbit pthread dl)

set(EXPERIMENT_PERFORMANCE_INTERSECT_SINGLE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_intersect_single.cpp
        )
add_executable(ex_performance_intersect_single ${EXPERIMENT_PERFORMANCE_INTERSECT_SINGLE_SOURCE_FILES})
target_link_libraries(ex_performance_intersect_single fastbit pthread dl)

set(EXPERIMENT_PERFORMANCE_DECOMPRESS_SINGLE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/common.hpp
        experiments/performance/main_performance_decompress_single.cpp
        )
add_executable(ex_performance_decompress_single ${EXPERIMENT_PERFORMANCE_DECOMPRESS_SINGLE_SOURCE_FILES})
target_link_libraries(ex_performance_decompress_single fastbit pthread dl)

# Index compression
set(EXPERIMENT_INDEX_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/common.hpp
        experiments/index/main_index_compression.cpp
        )
add_executable(ex_index_compression ${EXPERIMENT_INDEX_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_index_compression fastbit pthread dl)



set(EXPERIMENT_PERFORMANCE_SKYLINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_performance_skyline.cpp
        )
add_executable(ex_performance_skyline ${EXPERIMENT_PERFORMANCE_SKYLINE_SOURCE_FILES})
target_link_libraries(ex_performance_skyline fastbit pthread dl)

set(EXPERIMENT_PERFORMANCE_PROFILING_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_performance_profiling.cpp
        )
add_executable(ex_performance_profiling ${EXPERIMENT_PERFORMANCE_PROFILING_SOURCE_FILES})
target_link_libraries(ex_performance_profiling fastbit pthread dl)

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

set(EXPERIMENT_PERFORMANCE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/performance/main_performance.cpp
        )
add_executable(ex_performance ${EXPERIMENT_PERFORMANCE_SOURCE_FILES})
target_link_libraries(ex_performance fastbit)


set(EXPERIMENT_LOSSY_COMPRESSION_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/compression/main_lossycompression.cpp
        )
add_executable(ex_lossycompression ${EXPERIMENT_LOSSY_COMPRESSION_SOURCE_FILES})
target_link_libraries(ex_lossycompression fastbit pthread dl)

set(TOOL_COMPRESS_ALL_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/tool/compress_all.cpp
        )
add_executable(tool_compress_all ${TOOL_COMPRESS_ALL_SOURCE_FILES})
target_link_libraries(tool_compress_all fastbit pthread dl)

set(TOOL_HISTO_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/tool/histo.cpp
        )
add_executable(tool_histo ${TOOL_HISTO_SOURCE_FILES})
target_link_libraries(tool_histo fastbit pthread dl)

set(TOOL_COMPRESS_FILE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/tool/compress_file.cpp
        )
add_executable(compress_file ${TOOL_COMPRESS_FILE_SOURCE_FILES})
target_link_libraries(compress_file fastbit pthread dl)

set(REAL_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/real/real.cpp
        )
add_executable(real ${REAL_SOURCE_FILES})
target_link_libraries(real fastbit pthread dl boost_system boost_filesystem jemalloc)

set(REAL_LOSSY_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/real/real_lossy.cpp
        )
add_executable(real_lossy ${REAL_LOSSY_SOURCE_FILES})
target_link_libraries(real_lossy fastbit pthread dl boost_system boost_filesystem jemalloc)

set(REAL_LOSSY_EXAMINE_SOURCE_FILES
        ${SOURCE_FILES}
        ${BENCHMARK_SOURCE_FILES}
        experiments/real/real_lossy_examine.cpp
        )
add_executable(real_lossy_examine ${REAL_LOSSY_EXAMINE_SOURCE_FILES})
target_link_libraries(real_lossy_examine fastbit pthread dl boost_system boost_filesystem jemalloc)
##===----------------------------------------------------------------------===##
