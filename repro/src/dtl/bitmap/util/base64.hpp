#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap.hpp>
#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>

#include <iterator>
#include <string>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
//static std::vector<$u8>
//decode64(const std::string& val) {
//  using namespace boost::archive::iterators;
//  using It = transform_width<binary_from_base64<std::string::const_iterator>, 8, 6>;
//  return std::vector<$u8>(It(std::begin(val)), It(std::end(val)));
//}
//===----------------------------------------------------------------------===//
static std::vector<$u8>
decode64(const dtl::data_view<u8>& val) {
  using namespace boost::archive::iterators;
  using It = transform_width<
      binary_from_base64<typename std::vector<$u8>::const_iterator>, 8, 6>;
  return std::vector<$u8>(It(std::begin(val)), It(std::end(val)));
}
//===----------------------------------------------------------------------===//
static std::string
encode64(const dtl::data_view<u8>& val) {
  using namespace boost::archive::iterators;
  using It = base64_from_binary<
      transform_width<typename std::vector<$u8>::const_iterator, 6, 8>>;
  const auto b = std::begin(val);
  const auto e = std::end(val);
  return std::string(It(b), It(e));
}
//===----------------------------------------------------------------------===//
static std::string
base64_encode_bitmap(const dtl::bitmap& b) {
  u32 bitmap_length = b.size();
  std::vector<$u32> bitmap_data(b.num_blocks() + 1, 0u);
  bitmap_data[0] = bitmap_length;
  boost::to_block_range(b, bitmap_data.begin() + 1);
  const dtl::data_view<u8> bitmap_view {
    reinterpret_cast<u8*>(&bitmap_data[0]),
    reinterpret_cast<u8*>(&bitmap_data[0] + bitmap_data.size())
  };
  return encode64(bitmap_view);
}
//===----------------------------------------------------------------------===//
//static dtl::bitmap
//base64_decode_bitmap(const std::string& s) {
//  const auto decoded_bytes = decode64(s);
//  assert(decoded_bytes.size() % 4 == 0);
//  const dtl::data_view<u32> view {
//      reinterpret_cast<u32*>(&decoded_bytes[0] + 4),
//      reinterpret_cast<u32*>(&decoded_bytes[0] + decoded_bytes.size())
//  };
//
//  dtl::bitmap bm(view.begin(), view.end());
//  // Adjust the size.
//  u32 bitmap_length = *reinterpret_cast<u32*>(&decoded_bytes[0]);
//  bm.resize(bitmap_length);
//  return bm;
//}
//===----------------------------------------------------------------------===//
static dtl::bitmap
base64_decode_bitmap(const dtl::data_view<u8>& d) {
  const auto decoded_bytes = decode64(d);
  assert(decoded_bytes.size() % 4 == 0);
  const dtl::data_view<u32> view {
    reinterpret_cast<u32*>(&decoded_bytes[0] + 4),
    reinterpret_cast<u32*>(&decoded_bytes[0] + decoded_bytes.size())
  };
  dtl::bitmap bm(view.begin(), view.end());
  // Adjust the size.
  u32 bitmap_length = *reinterpret_cast<u32*>(&decoded_bytes[0]);
  bm.resize(bitmap_length);
  return bm;
}
//===----------------------------------------------------------------------===//
} // namespace dtl
//===----------------------------------------------------------------------===//
