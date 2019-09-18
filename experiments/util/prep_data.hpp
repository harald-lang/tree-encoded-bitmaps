#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_db.hpp"
#include "params.hpp"

#include <cstddef>
#include <vector>
//===----------------------------------------------------------------------===//
/// Ensures that 'cnt' random bitmaps for each configuration are stored in the
/// database. The function first consults the database, and generates new
/// bitmaps only if they are not present.
/// Note that invalid configurations are silently skipped.
void
prep_data(
    std::vector<params_markov>& params,
    std::size_t cnt, // the number of bitmaps to generate for each setting (typically corresponds to the number of independent runs)
    bitmap_db& db // the database, where the bitmaps are stored
);
void
prep_data(
    std::vector<params_uniform>& params,
    std::size_t cnt, // the number of bitmaps to generate for each setting (typically corresponds to the number of independent runs)
    bitmap_db& db // the database, where the bitmaps are stored
);
//===----------------------------------------------------------------------===//
