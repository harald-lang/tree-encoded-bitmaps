#pragma once

#include <dtl/adept.hpp>
#include <dtl/index.hpp>
#include <dtl/mem_info.hpp>
#include "histograms.hpp"
#include <boost/dynamic_bitset.hpp>

namespace col_imp{

#define MAXOFFSET 24

template<typename T, u32 BINS_CNT, bool EQ_W, bool block_wise, bool dict_comp = true>
class column_imprints{
  static_assert(BINS_CNT > 0 && !(BINS_CNT & (BINS_CNT-1)),
                "Template parameter 'BINS_CNT' must be a power of two.");

  struct dict_entry{
    $u64 cnt:MAXOFFSET; // number of imprints of this entry
    $u64 rpt:1; // repeat flag
    $u64 flgs:8 * sizeof(int) - MAXOFFSET-1; // for future use, e.g. micro swaps
  };

  dict_entry *dict;
  $u64 dict_num_entr;
  $u64 *imprints;
  $u64 imprints_num;
  static constexpr u32 gran = 64/sizeof(T);
  histo::histogram<T,BINS_CNT,EQ_W,block_wise> hist;
  size_t tpl_num;

public:
  struct col_imp_builder {
    column_imprints& ref;

    inline void
    create_imprints(const T *const values, const size_t n, const histo::histogram<T,BINS_CNT,EQ_W,block_wise>& h){
      $u64 i, mask, prev_mask;
      $u32 new_entry; // indicates if an dictionary entry is full
      ref.hist = h;
      ref.tpl_num = n;

      $u32 fits = 64 / ref.hist.bins_num_max; //how many mask vectors fit in the fingerprint
      ref.dict_num_entr = 0;
      ref.imprints_num = 0;
      prev_mask = mask = 0;

      $u64 tpl_per_impr = gran > n ? n : gran;
      new_entry = tpl_per_impr - 1; // tpl_per_impr is always power of 2
      $u64 num_max_impr = std::ceil(n / static_cast<$f64>(tpl_per_impr)); // max needed number of imprints

      ref.dict = (dict_entry*) malloc(sizeof(dict_entry) * num_max_impr);
      ref.imprints = ($u64*) malloc(sizeof($u64) * (std::ceil(static_cast<$f64>(num_max_impr) / fits)));

      //init imprints with 0s
      for (i = 0, ref.imprints_num = ((num_max_impr / fits) +1); i < ref.imprints_num; i++)
        ref.imprints[i] = 0;

      ref.imprints_num = 0;

      // start building the imprints
      for (i = 0; i < n; i++) {

        mask |= (1 << ref.hist.get_slot(values[i]));
        /*
        cout << "Slot: " << ref.hist.get_slot(values[i]) << endl;
        cout << "add new value: " << values[i] << ", bin: " << ref.hist.get_slot(values[i]) << " : "
             << std::bitset<BINS_CNT>(mask) << " shifted by: "
             << ((ref.imprints_num % fits) * ref.hist.bins_num_max) << " bits" << endl;
        cout << new_entry << " " << i << endl;
        */
        if (!(i+1 & new_entry) && i > 0) {  //aggregates all tuples of one cacheline
          //cout << "-----" << endl;
          if (prev_mask == mask && ref.dict[ref.dict_num_entr - 1].cnt < ((1 << MAXOFFSET) - 1)) { // compress list
            if (ref.dict[ref.dict_num_entr - 1].rpt == 0) {
              if (ref.dict[ref.dict_num_entr - 1].cnt > 1) {
                ref.dict[ref.dict_num_entr - 1].cnt--; // reduce previous by 1
                ref.dict_num_entr++;
                ref.dict[ref.dict_num_entr - 1].cnt = 1; // the new is a repeat
              }
              ref.dict[ref.dict_num_entr - 1].rpt = 1;
            }
            // same mask as before
            ref.dict[ref.dict_num_entr - 1].cnt++;
          } else {
            // new mask
            prev_mask = mask;
            ref.hist.bins_num_max == 64 ?
            (ref.imprints[ref.imprints_num] = mask) :
            (ref.imprints[ref.imprints_num / fits] |= mask << ((ref.imprints_num % fits) * ref.hist.bins_num_max));
            ref.imprints_num++;

            if(ref.dict_num_entr > 0 && ref.dict[ref.dict_num_entr - 1].rpt == 0 &&
               ref.dict[ref.dict_num_entr - 1].cnt < ((1 << MAXOFFSET) - 1)) {
              ref.dict[ref.dict_num_entr - 1].cnt++;
            } else {
              ref.dict[ref.dict_num_entr].cnt = 1;
              ref.dict[ref.dict_num_entr].rpt = 0;
              ref.dict_num_entr++;
            }
          }
          mask = 0;
          /*
          cout << "Dic_entry: " << ref.dict_num_entr << endl;
          cout << "      Cnt: " << ref.dict[ref.dict_num_entr-1].cnt << endl;
          cout << "   Repeat: " << ref.dict[ref.dict_num_entr-1].rpt << endl;
          */
        }
      }

      if((i+1 & new_entry)){
        // last mask, needed if (#tuples % gran) != 0
        if (prev_mask == mask && ref.dict_num_entr > 0 && ref.dict[ref.dict_num_entr - 1].cnt < (1 << MAXOFFSET)) {
          if (ref.dict[ref.dict_num_entr - 1].rpt == 0) {
            if (ref.dict[ref.dict_num_entr - 1].cnt == 1) // only 1 on previous
              ref.dict[ref.dict_num_entr - 1].rpt = 1;
            else {
              ref.dict[ref.dict_num_entr - 1].cnt--; // reduce previous by 1
              ref.dict[ref.dict_num_entr].cnt = 1; // the new is a repeat
              ref.dict[ref.dict_num_entr].rpt = 1;
              ref.hist.bins_num_max == 64 ? (ref.imprints[ref.imprints_num] = mask) :
              (ref.imprints[ref.imprints_num / fits] |= mask << ((ref.imprints_num % fits) * ref.hist.bins_num_max));
              ref.imprints_num++;
              ref.dict_num_entr++;
            }
          }
          // same mask as before
          ref.dict[ref.dict_num_entr - 1].cnt++;
        } else {
          ref.hist.bins_num_max == 64 ? (ref.imprints[ref.imprints_num] = mask) :
          (ref.imprints[ref.imprints_num / fits] |= mask << ((ref.imprints_num % fits) * ref.hist.bins_num_max));
          ref.imprints_num++;

          if (ref.dict_num_entr > 0 && ref.dict[ref.dict_num_entr - 1].rpt == 0
                                    && ref.dict[ref.dict_num_entr - 1].cnt < ((1 << MAXOFFSET) - 1)) {
            ref.dict[ref.dict_num_entr - 1].cnt++;
          } else {
            ref.dict[ref.dict_num_entr].cnt = 1;
            ref.dict[ref.dict_num_entr].rpt = 0;
            ref.dict_num_entr++;
          }
        }
      }
    }

    inline void
    create_imprints_no_dic_compr(const T *const values, const size_t n, const histo::histogram<T,BINS_CNT,EQ_W,block_wise>& h){
      $u64 i, mask, prev_mask;
      $u32 new_entry; // indicates if an dictionary entry is full
      ref.hist = h;
      ref.tpl_num = n;

      $u32 fits = 64 / ref.hist.bins_num_max; //how many mask vectors fit in the fingerprint
      ref.dict_num_entr = 0;
      ref.imprints_num = 0;
      mask = 0;

      $u64 tpl_per_impr = gran > n ? n : gran;
      new_entry = tpl_per_impr - 1; // tpl_per_impr is always power of 2
      $u64 num_max_impr = std::ceil(n / static_cast<$f64>(tpl_per_impr)); // max needed number of imprints

      ref.dict = (dict_entry*) malloc(sizeof(dict_entry) * num_max_impr);
      ref.imprints = ($u64*) malloc(sizeof($u64) * (std::ceil(static_cast<$f64>(num_max_impr) / fits)));

      //init imprints with 0s
      for (i = 0, ref.imprints_num = ((num_max_impr / fits) +1); i < ref.imprints_num; i++)
        ref.imprints[i] = 0;

      ref.imprints_num = 0;

      // start building the imprints
      for (i = 0; i < n; i++) {
        mask |= (1 << ref.hist.get_slot(values[i]));
/*
        cout << "Slot: " << ref.hist.get_slot(values[i]) << endl;
        cout << "add new value: " << values[i] << ", bin: " << ref.hist.get_slot(values[i]) << " : "
             << std::bitset<BINS_CNT>(mask) << " shifted by: "
             << ((ref.imprints_num % fits) * ref.hist.bins_num_max) << " bits" << endl;
        cout << new_entry << " " << i << endl;
*/
        if (!(i+1 & new_entry) && i > 0) {  //aggregates all tuples of one cacheline
          //cout << "-----" << endl;

          // new mask
          ref.hist.bins_num_max == 64 ?
          (ref.imprints[ref.imprints_num] = mask) :
          (ref.imprints[ref.imprints_num / fits] |= mask << ((ref.imprints_num % fits) * ref.hist.bins_num_max));
          ref.imprints_num++;

          if(ref.dict_num_entr > 0 && ref.dict[ref.dict_num_entr - 1].cnt < ((1 << MAXOFFSET) - 1)) {
            ref.dict[ref.dict_num_entr - 1].cnt++;
          } else {
            ref.dict[ref.dict_num_entr].cnt = 1;
            ref.dict[ref.dict_num_entr].rpt = 0;
            ref.dict_num_entr++;
          }

          mask = 0;
          /*
          cout << "Dic_entry: " << ref.dict_num_entr << endl;
          cout << "      Cnt: " << ref.dict[ref.dict_num_entr-1].cnt << endl;
          cout << "   Repeat: " << ref.dict[ref.dict_num_entr-1].rpt << endl;
          */
        }
      }

      if((i+1 & new_entry)){
        ref.hist.bins_num_max == 64 ? (ref.imprints[ref.imprints_num] = mask) :
        (ref.imprints[ref.imprints_num / fits] |= mask << ((ref.imprints_num % fits) * ref.hist.bins_num_max));
        ref.imprints_num++;

        if (ref.dict_num_entr > 0 && ref.dict[ref.dict_num_entr - 1].cnt < ((1 << MAXOFFSET) - 1)) {
          ref.dict[ref.dict_num_entr - 1].cnt++;
        } else {
          ref.dict[ref.dict_num_entr].cnt = 1;
          ref.dict[ref.dict_num_entr].rpt = 0;
          ref.dict_num_entr++;
        }
      }
    }

    inline void
    operator()(const T *const values, const size_t n){

      histo::samples<T> s = histo::samples<T>(values, n);
      histo::histogram<T,BINS_CNT,EQ_W,block_wise> h(s);

    }

    inline void
    operator()(const T *const values, const size_t n, const histo::histogram<T,BINS_CNT,EQ_W,block_wise>& h){

      if(dict_comp)
        create_imprints(values, n, h);
      else
        create_imprints_no_dic_compr(values,n,h);
    }

    inline void done(){
      std::free(ref.dict);
      std::free(ref.imprints);
    }

    void print(){
      $u64 impr = 0;
      $u64 count = 0;
      T* imprints = reinterpret_cast<T*>(ref.imprints);

      auto get_impr = [&](const auto num){
        if(BINS_CNT == 64)
          return std::bitset<BINS_CNT>(ref.imprints[num]);
        else{
          $u32 fits = 64 / ref.hist.bins_num_max;
          $u32 shift = (num % fits) * ref.hist.bins_num_max;
          return std::bitset<BINS_CNT>(ref.imprints[num / fits] >> shift);
        }
      };

      for(auto i = 0; i < ref.dict_num_entr; i++){
        std::cout << "-----" << std::endl;
        std::cout << "Dic_entry: " << i << std::endl;
        std::cout << "      Cnt: " << ref.dict[i].cnt << std::endl;
        std::cout << "   Repeat: " << ref.dict[i].rpt << std::endl;
        std::cout << "-" << std::endl;

        if(ref.dict[i].rpt) {
          for(auto j = 0; j < ref.dict[i].cnt; j++)
            std::cout << count++ << ": " << get_impr(impr).to_string() << std::endl;

          impr++;
        } else
          for(auto j = 0; j < ref.dict[i].cnt; j++)
            std::cout << count++ << ": " << get_impr(impr++).to_string() << std::endl;

        std::cout << "-" << std::endl;
      }
    }


    inline dtl::mem_info memory_footprint(){
      // TODO memory_consumption
      auto bytes_ptr = sizeof(dict) + sizeof(imprints);
      auto bytes_dict = sizeof(dict_entry) * ref.dict_num_entr;
      auto bytes_impr = sizeof(*ref.imprints) * ref.imprints_num;
      auto bytes_histo = ref.hist.memory_footprint();
      auto bytes_local_var = sizeof(gran) + sizeof(tpl_num);

      std::string idx;
      if(dict_comp)
        idx = "Column Imprints";
      else
        idx = "Column Imprints without dictionary compression";

      struct dtl::mem_info memory(idx, BINS_CNT, bytes_ptr, bytes_impr + bytes_dict, bytes_local_var, bytes_histo);

      return memory;
    }

  };

  inline col_imp_builder
  builder() {
    // return a builder instance
    return col_imp_builder { *this };
  }

  // c'tor
  column_imprints() noexcept {
    // initialize the variables

  }

  template<typename I>
  inline boost::dynamic_bitset<> lookup_helper(boost::dynamic_bitset<>& r, $u8 b, u8 e){
    std::cout << "use helper" << std::endl;
    std::cout << "mask: set bits [" << std::to_string(b) << ";" << std::to_string(e) << "]" << std::endl;

    I* imp_ptr = reinterpret_cast<I*>(imprints);
    $u64 imp_pos = 0;
    $u64 tpl_pos = 0;
    I mask = 0;

    for(auto i = b; i <= e; i++)
      mask |= (1 << i);

    std::cout << "mask: " << std::bitset<BINS_CNT>(mask) << std::endl;

    auto set_cl = [&](){
      for(auto i = 0; i < gran; i++)
        r.set(tpl_pos+i, true);
      tpl_pos += gran;
    };

    for(auto i = 0; i < dict_num_entr; i++){
      if(dict[i].rpt){
        if((imp_ptr[imp_pos++] & mask) > 0)
          for(auto c = 0; c < dict[i].cnt; c++)
            set_cl();
        else
          tpl_pos += dict[i].cnt*gran;

      }else{
        for(auto c = 0; c < dict[i].cnt; c++) {
          if((imp_ptr[imp_pos++] & mask) > 0)
            set_cl();
          else
            tpl_pos += gran;
        }
      }
    }

    return r;
  }

  inline boost::dynamic_bitset<>
  lookup(const dtl::predicate& p){
    using value_t = typename std::remove_cv<T>::type;

    boost::dynamic_bitset<> r(tpl_num);
    value_t value = *reinterpret_cast<value_t*>(p.value_ptr);
    value_t second_value; // in case of between predicates

    if(block_wise){
      if(value < hist.get_min()){
        switch (p.comparison_operator) {
          case dtl::op::EQ:
          case dtl::op::LE:
          case dtl::op::LT:
            return r;
          case dtl::op::GT:
          case dtl::op::GE:
            r.set();
            return r;
        }
      }

      if(value > hist.get_max()){
        switch (p.comparison_operator){
          case dtl::op::EQ:
          case dtl::op::GT:
          case dtl::op::GE:
            return r;
          case dtl::op::LE:
          case dtl::op::LT:
            r.set();
            return r;
        }
      }
    }

    $u8 b = 0;
    $u8 e = 0;

    $u32 bin_value = hist.get_slot(value);

    switch(p.comparison_operator){
      case dtl::op::EQ:
        b = e = bin_value;
        break;

      case dtl::op::LT:
        if(value == hist.get_min())
          return r;
        e = hist.get_slot(histo::next_smaller_value<T>(value));
        break;

      case dtl::op::LE:
        e = hist.get_slot(value);
        break;

      case dtl::op::GT:
        if(value == hist.get_max())
          return r;
        b = hist.get_slot(histo::next_greater_value<T>(value));
        e = BINS_CNT-1;
        break;

      case dtl::op::GE:
        b = hist.get_slot(value);
        e = BINS_CNT-1;
        break;

      case dtl::op::BETWEEN:
        if(value > hist.get_max())
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < hist.get_min())
          return r;

        if(value < hist.get_min())
          b = 0;
        else
          b = hist.get_slot(value);

        e = hist.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_LO: // (x,y]
        if(value >= hist.get_max()) // the searched value is the largest possible value => no larger values
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < hist.get_min()) // upper bound smaller than smallest value
          return r;

        if(value < hist.get_min()) // first value smaller than minimum
          b = 0;
        else
          b = hist.get_slot(histo::next_greater_value(value)); // exclude x

        e = hist.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_RO: // [x,y)
        if(value > hist.get_max())
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= hist.get_min()) // upper bound smaller than smallest value
          return r;

        if(value < hist.get_min()) // first value smaller than minimum
          b = 0;
        else
          b = hist.get_slot(value);

        e = hist.get_slot(histo::next_smaller_value(second_value)); // exclude y
        break;

      case dtl::op::BETWEEN_O: // (x,y)
        if(value >= hist.get_max())
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= hist.get_min())
          return r;

        if(value < hist.get_min())
          b = 0;
        else
          b = hist.get_slot(histo::next_greater_value(value)); // exclude x -> [x+1;y)

        e = hist.get_slot(histo::next_smaller_value(second_value)); // exclude y -> [x+1;y-1]
        break;
    }

    switch(BINS_CNT){
      case 8:
        return lookup_helper<$u8>(r, b, e);
      case 16:
        return lookup_helper<$u16>(r, b, e);
      case 32:
        return lookup_helper<$u32>(r, b, e);
      case 64:
        return lookup_helper<$u64>(r, b, e);
    }
  }
};

} // namespace col_imp



