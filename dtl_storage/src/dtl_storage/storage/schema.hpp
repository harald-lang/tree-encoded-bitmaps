#pragma once

#ifndef _DTL_STORAGE_INCLUDED
#error "Never use <dtl/storage/schema.hpp> directly; include <dtl/storage.hpp> instead."
#endif

#include <string>
#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/storage/types.hpp>


namespace dtl {


struct attr {
  std::string name;
  dtl::rtt type;
};

using schema = std::vector<attr>;


} // namespace dtl