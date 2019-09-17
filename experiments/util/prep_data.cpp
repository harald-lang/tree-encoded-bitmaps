#include "prep_data.hpp"

#include <dtl/dtl.hpp>

#include "gen.hpp"
#include "params.hpp"
//===----------------------------------------------------------------------===//
void
prep_data(
    std::vector<params_markov>& params,
    const std::size_t cnt,
    bitmap_db& db
) {

  // Prepare the random bitmaps.
  std::cout << "Preparing the data set." << std::endl;
  std::vector<params_markov> missing_bitmaps;

  for (auto& p : params) {
    const auto n = p.n;
    const auto f = p.clustering_factor;
    const auto d = p.density;
    if (!markov_parameters_are_valid(n, f, d)) continue;

    std::cout << n << "," << d << "," << f << std::endl;
    auto ids = db.find_bitmaps(n, f, d);
    if (ids.size() < cnt) {
      params_markov c;
      c.n = n;
      c.clustering_factor = f;
      c.density = d;

      for (std::size_t i = ids.size(); i < cnt; ++i) {
        missing_bitmaps.push_back(c);
      }
    }
  }

  if (!missing_bitmaps.empty()) {
    std::cout << "Generating " << missing_bitmaps.size()
        << " random bitmaps." << std::endl;

    const std::size_t failure_cntr = gen(missing_bitmaps, db);

    if (failure_cntr > 0) {
      std::cerr << "Failed to generate all required bitmaps. "
          << failure_cntr << " bitmaps are still missing."
          << std::endl;
    }
  }

  std::size_t pass = 2;
  while (true) {
    std::cout << "Preparing the data set. (pass " << pass << ")" << std::endl;
    std::vector<params_markov> incomplete_bitmaps;
    for (auto& p : params) {
      const auto n = p.n;
      const auto f = p.clustering_factor;
      const auto d = p.density;

      if (!markov_parameters_are_valid(n, f, d)) continue;

      auto ids = db.find_bitmaps(n, f, d);
      if (ids.size() > 0 && ids.size() < cnt) {
        params_markov c;
        c.n = n;
        c.clustering_factor = f;
        c.density = d;
        for (std::size_t i = ids.size(); i < cnt; ++i) {
          incomplete_bitmaps.push_back(c);
        }
      }

    }
    std::cout << incomplete_bitmaps.size() << " remaining." << std::endl;
    if (!incomplete_bitmaps.empty()) {
      std::cout << "Generating " << incomplete_bitmaps.size()
          << " random bitmaps. (pass " << pass << ")" << std::endl;

      const std::size_t failure_cntr = gen(incomplete_bitmaps, db);

      if (failure_cntr == 0) {
        break;
      }

    }
    else {
      break;
    }
    pass++;
  }
  std::cerr << "Done generating random bitmaps after "
      << pass << " passes." << std::endl;
}
//===----------------------------------------------------------------------===//
void
prep_data(
    std::vector<params_uniform>& params,
    const std::size_t cnt,
    bitmap_db& db
) {

  // Prepare the random bitmaps.
  std::cout << "Preparing the data set." << std::endl;
  std::vector<params_uniform> missing_bitmaps;

  for (auto& p : params) {
    const auto n = p.n;
    const auto d = p.density;
    if (n == 0 || d < 0 || d > 1.0) continue;

    auto ids = db.find_bitmaps(n, d);
    if (ids.size() < cnt) {
      params_uniform c;
      c.n = n;
      c.density = d;

      for (std::size_t i = ids.size(); i < cnt; ++i) {
        missing_bitmaps.push_back(c);
      }
    }
  }

  if (!missing_bitmaps.empty()) {
    std::cout << "Generating " << missing_bitmaps.size()
        << " random bitmaps." << std::endl;

    const std::size_t failure_cntr = gen(missing_bitmaps, db);

    if (failure_cntr > 0) {
      std::cerr << "Failed to generate all required bitmaps. "
          << failure_cntr << " bitmaps are still missing."
          << std::endl;
    }
  }

  std::size_t pass = 2;
  while (true) {
    std::cout << "Preparing the data set. (pass " << pass << ")" << std::endl;
    std::vector<params_uniform> incomplete_bitmaps;
    for (auto& p : params) {
      const auto n = p.n;
      const auto d = p.density;
      if (n == 0 || d < 0 || d > 1.0) continue;

      auto ids = db.find_bitmaps(n, d);
      if (ids.size() > 0 && ids.size() < cnt) {
        params_uniform c;
        c.n = n;
        c.density = d;
        for (std::size_t i = ids.size(); i < cnt; ++i) {
          incomplete_bitmaps.push_back(c);
        }
      }

    }
    std::cout << incomplete_bitmaps.size() << " remaining." << std::endl;
    if (!incomplete_bitmaps.empty()) {
      std::cout << "Generating " << incomplete_bitmaps.size()
          << " random bitmaps. (pass " << pass << ")" << std::endl;

      const std::size_t failure_cntr = gen(incomplete_bitmaps, db);

      if (failure_cntr == 0) {
        break;
      }

    }
    else {
      break;
    }
    pass++;
  }
  std::cerr << "Done generating random bitmaps after "
      << pass << " passes." << std::endl;
}
//===----------------------------------------------------------------------===//
