#pragma once

#include <boost/tokenizer.hpp>
#include <dtl/storage.hpp>

namespace fp{

  std::vector<std::string>
  parse_csv_line(const std::string& csv_input,
                 const std::string separator) {

    using namespace boost;
    boost::char_separator<char> sep(separator.c_str(), "", boost::keep_empty_tokens);
    typedef tokenizer<char_separator<char>> Tokenizer;

    std::vector<std::string> tokens;

    Tokenizer tok(csv_input, sep);
    tokens.assign(tok.begin(), tok.end());
    return tokens;
  }

// parse CSV input and create column blocks
  template<u64 block_size>
  std::vector<dtl::column_block_base<block_size>*>
  parse_input_tuples(const dtl::schema s, const std::vector<std::string>& csv) {
    if (csv.size() != block_size) {
      std::cout << "!" << std::flush;
    }
    std::string null_indicator("NA");

    // create a column block for each attribute
    std::vector<dtl::column_block_base<block_size>*> blocks;
    for (auto& attr : s) {
      blocks.push_back(dtl::make_column_block<block_size>(attr.type));
    }

    // parse input
    for (auto& csv_line : csv) {
      auto values = parse_csv_line(csv_line, ",");

      // insert attribute values
      for ($u64 i = 0; i < s.size(); i++) {
        auto block_ptr = blocks[i];
        const dtl::rtt type = s[i].type;
        try {
          dtl::column_block_insert<block_size>(block_ptr, type, values[i], null_indicator);
        }
        catch(...) {
          std::cerr << "Failed to parse value '" << values[i] << "' for attribute " << i << "." << std::endl;
        }
      }
    }

    return blocks;
  }

};

