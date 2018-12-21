#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/env.hpp>
#include <dtl/thread.hpp>

/// The number of independent runs.
static constexpr std::size_t REPEAT = 10;

using tree_t = dtl::binary_tree_structure;

struct setting_t {
  $u64 n = 0;
  $f64 d = 0.0;
  $f64 f = 0.0;
};

void
run_verbose(const setting_t& setting) {
  u64 n = setting.n;
  f64 f = setting.f;
  f64 d = setting.d;

  const auto plain_bitmap = dtl::gen_random_bitmap(n, f, d);
  std::cout << "plain bitmap:        " << plain_bitmap << std::endl;

  const std::vector<$f64> fprs {0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0};
//  const std::vector<$f64> fprs {0.1};
  for (auto fpr : fprs) {
    std::cout << "desired fpr: " << fpr << std::endl;
    const auto compressed_bitmap = dtl::dynamic_tree_mask_lo(plain_bitmap, fpr);
    const auto decompressed_bitmap = compressed_bitmap.to_bitset();
    const auto actual_false_positive_cnt =
        (decompressed_bitmap ^ plain_bitmap).count();

    const auto tp_cnt = plain_bitmap.count();
    const auto tn_cnt = n - tp_cnt;
    u64 max_fp_cnt = static_cast<u64>(tn_cnt * fpr);
    std::cout << "max fp cnt:  " << max_fp_cnt << std::endl;
    const auto actual_fpr = (actual_false_positive_cnt * 1.0) / n;
    std::cout << "actual fpr:  "
              << actual_false_positive_cnt << " / " << n << " = "
              << actual_fpr << std::endl;

    std::cout << "compressed bitmap:   " << compressed_bitmap << std::endl;
    std::cout << "decompressed bitmap: " << decompressed_bitmap << std::endl;
    std::cout << "compression ratio:   "
        << ((compressed_bitmap.serialized_size_in_byte() * 1.0) / ((n + 7) / 8))
        << std::endl;
    {
      // Validation code.
      for (std::size_t i = 0; i < n; ++i) {
        if (plain_bitmap[i] && !decompressed_bitmap[i]) {
          std::cout << "Validation failed." << std::endl;
          std::exit(1);
        }
      }
    }
    std::cout << std::endl;

  }
}

void
run(const setting_t& setting) {
  u64 n = setting.n;
  f64 f = setting.f;
  f64 d = setting.d;

  const auto plain_bitmap = dtl::gen_random_bitmap(n, f, d);

  const std::vector<$f64> fprs {0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0};
  for (auto fpr : fprs) {
    const auto compressed_bitmap = dtl::dynamic_tree_mask_lo(plain_bitmap, fpr);
    const auto decompressed_bitmap = compressed_bitmap.to_bitset();
    const auto actual_false_positive_cnt =
        (decompressed_bitmap ^ plain_bitmap).count();

    const auto tp_cnt = plain_bitmap.count();
    const auto tn_cnt = n - tp_cnt;
    const auto actual_fpr = (actual_false_positive_cnt * 1.0) / n;

    {
      // Validation code.
      if (actual_fpr > fpr) {
        std::cerr << "Validation failed. Actual FPR is higher than expected."
                  << std::endl;
        std::exit(1);
      }

      for (std::size_t i = 0; i < n; ++i) {
        if (plain_bitmap[i] && !decompressed_bitmap[i]) {
          std::cerr << "Validation failed. Found false negative results."
                    << std::endl;
          std::exit(1);
        }
      }
    }

    dtl::dynamic_roaring_bitmap compressed_with_roaring(plain_bitmap);
    dtl::dynamic_roaring_bitmap lossy_compressed_with_roaring(decompressed_bitmap);
    dtl::dynamic_wah32 compressed_with_wah(plain_bitmap);
    dtl::dynamic_wah32 lossy_compressed_with_wah(decompressed_bitmap);
    std::stringstream out;

    // Determine the actual bit density and the clustering factor.
    const auto actual_d = (plain_bitmap.count() * 1.0) / n;
    const auto actual_f =
        (plain_bitmap.count() * 1.0) / dtl::count_1fills(plain_bitmap);

    const auto compressed_bitmap_actual_d =
        (decompressed_bitmap.count() * 1.0) / n;
    const auto compressed_bitmap_actual_f =
        (decompressed_bitmap.count() * 1.0) / dtl::count_1fills(decompressed_bitmap);

    out << n
        << "," << d
        << "," << actual_d
        << "," << f
        << "," << actual_f
        << "," << fpr
        << "," << actual_fpr
        << "," << compressed_bitmap.structure_.size()
        << "," << compressed_bitmap.labels_.size()
        << "," << compressed_bitmap_actual_d
        << "," << compressed_bitmap_actual_f
        << "," << compressed_with_roaring.size_in_byte()
        << "," << lossy_compressed_with_roaring.size_in_byte()
        << "," << compressed_with_wah.size_in_byte()
        << "," << lossy_compressed_with_wah.size_in_byte()
        << std::endl;
    std::cout << out.str();
  }
}

$i32 main() {
  run_verbose(setting_t{.n = 1u << 10, .d = 0.1, .f = 16});
//  return 0;

//  const std::vector<$u64> n_values { 1u << 10, 1u << 12, 1u << 14, 1u << 16, 1u << 18, 1u << 20 };
  const std::vector<$u64> n_values { 1u << 20 };
  const std::vector<$f64> densities {0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1};

  std::vector<setting_t> settings;
  for(auto f = 1; f < 400; f += 1) {
    for (auto d : densities) {
      for (auto n : n_values) {
        settings.emplace_back(
            setting_t{.n = n, .d = d, .f = static_cast<f64>(f)});
      }
    }
  }
  std::random_shuffle(settings.begin(), settings.end());

  const auto config_cnt = settings.size();
  const std::size_t min_inc = 1;//32;
  const auto time_start = std::chrono::system_clock::now();
  std::atomic<std::size_t> cntr { 0 };
  auto thread_fn = [&](u32 thread_id) {
    while (true) {
      // Grab work.
      const auto inc = min_inc;
      const std::size_t config_idx_begin = cntr.fetch_add(inc);
      const std::size_t config_idx_end = std::min(config_idx_begin + inc, config_cnt);
      if (config_idx_begin >= config_cnt) break;

      for (std::size_t ci = config_idx_begin; ci < config_idx_end; ci++) {
        for (std::size_t r = 0; r < REPEAT; ++r){
          try {
            run(settings[ci]);
          }
          catch (...) {
            std::cerr << "Failed to construct bitmap." << std::endl;
          }
        }
      }

      if (thread_id == 0) {
        u64 i = cntr;
        u64 r = std::min(config_cnt, config_cnt - i);
        // Estimate time until completion.
        const auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = now - time_start;
        f64 avg_sec_per_config = elapsed_seconds.count() / i;
        u64 remaining_sec = avg_sec_per_config * r;
        u64 h = (remaining_sec / 3600);
        u64 m = (remaining_sec % 3600) / 60;
        std::stringstream str;
        str << "Progress of current run: [" << (i + 1) << "/" << config_cnt << "]";
        str << " - est. time until completion: " << h << "h " << m << "m" << std::endl;
        std::cerr << str.str();
      }
    }
  };
  dtl::run_in_parallel(thread_fn);


return 0;
//  auto t = dtl::dynamic_full_binary_tree(64);
//  for (std::size_t l = 0; l < t.height; ++l) {
//
//    $u64 from_node_idx = 0;
//    $u64 to_node_idx = 0;
//    for (std::size_t i = 0; i < l; ++i) {
//      from_node_idx = t.left_child_of(from_node_idx);
//      to_node_idx = t.right_child_of(to_node_idx);
//    }
//    assert(t.level_of(from_node_idx) == l);
//    assert(t.level_of(to_node_idx) == l);
//    std::cout << "level=" << l
//              << ": [" << from_node_idx << "," << to_node_idx << "]"
//              << ", #=" << (to_node_idx - from_node_idx + 1) << ""
//              << " - [" << (1ull << l) -1 << "," << (2*((1ull << l) - 1)) << "]"
//              << " - [" << (1ull << l) -1 << "," << (((1ull << (l+1)) - 2)) +1 << ")"
//              << std::endl;
//  }
}