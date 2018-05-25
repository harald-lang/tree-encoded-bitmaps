#pragma once

#include <dtl/dtl.hpp>
#include <algorithm>
#include <unordered_map>

namespace histo{

  template<typename T> T
  next_greater_value(T v){
    std::string v_type(typeid(v).name());

    if(v_type == "d")
      return nextafter(v, std::numeric_limits<T>::max());
    if(v_type == "f")
      return nextafter(v, std::numeric_limits<T>::max());
    if(v_type == "e")
      return nextafter(v, std::numeric_limits<T>::max());

    if(v == std::numeric_limits<T>::max())
      return v; // TODO wie loese ich das am besten?

    return ++v;
  }

  template<typename T> T
  next_smaller_value(T v){
    std::string v_type(typeid(v).name());

    if(v_type == "d")
      return nextafter(v, std::numeric_limits<T>::min());
    if(v_type == "f")
      return nextafter(v, std::numeric_limits<T>::min());
    if(v_type == "e")
      return nextafter(v, std::numeric_limits<T>::min());

    if(v == std::numeric_limits<T>::min())
      return v; // TODO wie loese ich das am besten?

    return --v;
  }

  u32 iter_eq_h = 5; // defines the #iterations for improving the eq_height-histogram

  template<typename T>
  struct value{
    T val;
    $u64 freq; // frequency of the value in the data

    value(const T& val, $u64 freq) : val(val), freq(freq){}

    void print(){
      std::cout << std::to_string(val) << " : " << freq << std::endl;
    }
  };

  template<typename T>
  struct samples{
    std::vector<value<T>> v;

    samples(const std::vector<T>& data){
      std::vector<T> dist_val = data;
      std::sort(dist_val.begin(), dist_val.end());
      auto last = std::unique(dist_val.begin(), dist_val.end());
      dist_val.erase(last, dist_val.end());

      for(auto e: dist_val)
        v.push_back(value<T>(e, std::count(data.begin(),data.end(), e)));
    }

    samples(const T* data, const size_t n){
      std::unordered_map<T, $u64> u;

      u.reserve(static_cast<$u64>(std::ceil(n/3.0)));

      for($u32 i = 0; i < n; i++)
        if (u.find(data[i]) != u.end())
          u[data[i]]++;
        else
          u[data[i]] = 1;

      for(const auto &e: u)
        v.push_back(value<T>(e.first, e.second));

      auto cmp = [&](const value<T> &t1, const value<T> &t2){
        return t1.val < t2.val ? true : false;
      };

      std::sort(v.begin(), v.end(), cmp);
    }

    samples(const std::vector<samples<T>>& s){
      std::unordered_map<T, $u64> u;
      $u64 res_size = 0;

      for(const auto &e: s)
        res_size += e.v.size();

      u.reserve(static_cast<$u64>(std::ceil(res_size/2.0)));

      for(const auto &e: s){
        for(auto d: e.v) {
          if (u.find(d.val) != u.end())
            u[d.val] += d.freq;
          else
            u[d.val] = d.freq;
        }
      }

      for(const auto &e: u)
        v.push_back(value<T>(e.first, e.second));

      auto cmp = [&](const value<T> &t1, const value<T> &t2){
        return t1.val < t2.val ? true : false;
      };

      std::sort(v.begin(), v.end(), cmp);
    }

    void print_min(){
      if(v.size() != 0){
        std::cout << "Min: ";
        v.front().print();
      }
      else
        std::cout << "Vector is empty! No min to print." << std::endl;
    }

    void print_max(){
      if(v.size() != 0) {
        std::cout << "Max: ";
        v.back().print();
      }else
        std::cout << "Vector is empty! No max to print." << std::endl;
    }

    void print(){
      if(v.size() != 0)
        for(auto e: v)
          e.print();
      else
        std::cout << "Vector is empty! Nothing to print." << std::endl;
    }
  };

  // TODO @Harald Bins als Array oder Vector -> vector waere dynamisch
  template<typename T, u32 BINS_CNT, bool EQ_W, bool block_wise>
  struct histogram{
    $u32 bins_used = 0;
    $u32 bins_num_max = BINS_CNT;
    T bin_min[BINS_CNT];
    T bin_max[BINS_CNT];

    histogram(){
    }

    histogram(const samples<T>& data){
      if(block_wise){
        if(data.v.size() < BINS_CNT -1)
          blockwise_direct(data);
        else{
          if(EQ_W)
            blockwise_eq_w(data);
          else
            blockwise_eq_h(data);
        }
      }else{
        if(data.v.size() < BINS_CNT -1)
          direct(data);
        else{
          if(EQ_W)
            eq_w(data);
          else
            eq_h(data);
        }
      }
    }

    histogram(const vector<T>& data) : histogram(samples<T>(data)){}

    histogram(const T* data, const size_t n) : histogram(samples<T>(data,n)){}

    histogram& operator=(histogram&& other){
      this->bins_used = other.bins_used;

      for(auto i = 0; i < BINS_CNT; i++){
        this->bin_min[i] = other.bin_min[i];
        this->bin_max[i] = other.bin_max[i];
      }
      return *this;
    }
    // TODO check
    histogram& operator=(const histogram& other){
      this->bins_used = other.bins_used;

      for(auto i = 0; i < BINS_CNT; i++){
        this->bin_min[i] = other.bin_min[i];
        this->bin_max[i] = other.bin_max[i];
      }
      return *this;
    }

    void eq_w(const samples<T>& data){
      u64 size_data = data.v.size();

      $u32 i; // for the current position
      bin_min[0] = std::numeric_limits<T>::min();
      bin_max[0] = data.v[0].val;
      bin_min[BINS_CNT - 1] = next_greater_value(data.v[size_data - 1].val);
      bin_max[BINS_CNT - 1] = std::numeric_limits<T>::max();

      $f64 y, y_step = static_cast<$f64>(size_data + 1) / static_cast<$f64>(BINS_CNT - 2);
      for(i = 1, y = y_step; y < size_data; y += y_step, i++) {
        bin_min[i] = bin_max[i - 1];
        bin_max[i] = data.v[static_cast<$u64>(y)].val;
      }
      // fill the last regular BIN [val; smallest_greater_value_than_max)
      bin_min[i] = bin_max[i-1];
      bin_max[i] = next_greater_value(data.v[size_data - 1].val);
      bins_used = i+2;
    }

    void eq_h(const samples<T>& data){
      $u64 height;
      $u64 k, cnt, j;
      $u64 usedBins;
      u64 size_data = data.v.size();

      height = size_data / (size_data < BINS_CNT - 2 ? size_data : (BINS_CNT - 2));

      bin_min[0] = std::numeric_limits<T>::min();
      bin_max[0] = data.v[0].val;
      bin_min[BINS_CNT - 1] = next_greater_value(data.v[size_data - 1].val);
      bin_max[BINS_CNT - 1] = std::numeric_limits<T>::max();

      $i64 normalizedCount[size_data];
      $i64 elements, normalized;
      for($u32 n = 0; n < iter_eq_h ; n++){
        elements = 0;
        normalized = 0;

        for($u64 i = 0; i < size_data; i++){
          if(data.v[i].freq > height) {
            normalizedCount[i] = 0;
            elements++;
          }else
            normalizedCount[i] = data.v[i].freq;
        }

        for($u64 i = 0; i < size_data; i++)
          normalized += normalizedCount[i];

        $u64 norm_h = ceil(normalized / (BINS_CNT - 2 - elements));
        if(height == norm_h)
          break;

        height = norm_h;
      }

      for(k = 1, j = 0; k < BINS_CNT - 1; k++) {
        bin_min[k] = bin_max[k - 1];
        cnt = data.v[j++].freq;
        while ((j < size_data) && (cnt < height)){
          if(cnt + data.v[j].freq >= 2*height)
            break;

          cnt += data.v[j++].freq;
        }
        bin_max[k] = data.v[j].val;
        if(j == size_data || k == BINS_CNT-2){
          bin_max[k] = next_greater_value(data.v[size_data-1].val);
          break;
        }
      }

      bins_used = k+2;

      for(k++; k < BINS_CNT - 1; k++) {
        bin_min[k] = next_greater_value(data.v[size_data-1].val);
        bin_max[k] = std::numeric_limits<T>::max();
      }
    }

    // create a histogram with only one value in each bin
    void direct(const samples<T>& data){
      u64 size_data = data.v.size();

      $u32 i = 1;

      if(size_data == 1) {
        bin_min[i] = data.v[0].val;
        bin_max[i] = next_greater_value(data.v[0].val);
      }else{
        for(i = 1; i < size_data; i++) {
          bin_min[i] = data.v[i - 1].val;
          bin_max[i] = data.v[i].val;
        }
        // create the bin with [max_sample; next_greater_val)
        bin_min[i] = data.v[i - 1].val;
        bin_max[i] = next_greater_value(bin_min[i]);
      }
      bins_used = i+1;

      for(int k = i+1; k < BINS_CNT; k++){
        bin_min[k] = bin_max[i];
        bin_max[k] = std::numeric_limits<T>::max();
      }

      //Fill the first BIN: [-inf;??)
      bin_min[0] = std::numeric_limits<T>::min();
      bin_max[0] = bin_min[1];
    }

    void blockwise_eq_w(const samples<T>& data){
      u64 size_data = data.v.size();

      bin_min[0] = data.v[0].val;
      bin_max[BINS_CNT-1] = next_greater_value(data.v[size_data-1].val);

      $f64 y, y_step = static_cast<$f64>(size_data) / static_cast<$f64>(BINS_CNT);
      $u32 i;
      for(i = 1, y = y_step; y < size_data; y += y_step, i++){
        bin_min[i] = data.v[static_cast<$u64>(y)].val;
        bin_max[i-1] = data.v[static_cast<$u64>(y)].val;
      }

      bins_used = i;
    }

    void blockwise_eq_h(const samples<T>& data){
      u64 size_data = data.v.size();

      $u64 h;
      $u64 k, cnt, j;
      $u64 tpl_cnt = 0;

      for(auto i = 0; i < size_data; i++)
        tpl_cnt += data.v[i].freq;

      h = tpl_cnt / BINS_CNT;

      $u64 norm_cnt[size_data];
      $u64 elem, norm;
      for($u16 n = 0; n < iter_eq_h ; n++){
        elem = 0;
        norm = 0;

        for($u64 i = 0; i < size_data; i++){
          if(data.v[i].freq >= h) {
            norm_cnt[i] = 0;
            elem++;
          }else
            norm_cnt[i] = data.v[i].freq;
        }
        for(auto i = 0; i < size_data; i++)
          norm += norm_cnt[i];

        $u64 norm_h = norm / (BINS_CNT - elem);

        if(h == norm_h)
          break;
        else
          h = norm_h;
      }

      for(k = 0, j = 0; k < BINS_CNT; k++) {
        bin_min[k] = k ? bin_max[k - 1] : data.v[0].val;
        cnt = data.v[j++].freq;

        while((j < size_data) && (cnt < h)){
          if(cnt + data.v[j].freq >= 2*h)
            break;

          cnt += data.v[j++].freq;
        }

        bin_max[k] = data.v[j].val;

        if(k == BINS_CNT-1 || j == size_data) {
          bin_max[k] = next_greater_value(data.v[size_data - 1].val);
          break;
        }
      }

      bins_used = k+1;

      for(k++; k < BINS_CNT; k++) {
        bin_min[k] = next_greater_value(data.v[size_data-1].val);
        bin_max[k] = std::numeric_limits<T>::max();
      }
    }

    // create a histogram with only one value in each bin
    void blockwise_direct(const samples<T>& data){
      $u16 k = 0;
      u64 size_data = data.v.size();

      if(size_data == 1) {
        bin_min[0] = data.v[0].val;
        bin_max[0] = next_greater_value(data.v[0].val);
      }else{
        for(k = 0; k < size_data - 1; k++) {
          bin_min[k] = data.v[k].val;
          bin_max[k] = data.v[k + 1].val;
        }
        bin_min[k] = data.v[size_data - 1].val;
        bin_max[k] = next_greater_value(data.v[size_data - 1].val);
      }
      bins_used = k+1;

      // Fill the remaining bins, needed for Column_Imprints, but not for PSMAs
      for(k++; k < BINS_CNT; k++) {
        bin_min[k] = next_greater_value(data.v[size_data - 1].val);
        bin_max[k] = next_greater_value(data.v[size_data - 1].val);
      }
    }

    void reset(){
      for(auto i = 0; i < BINS_CNT; i++){
        bin_min[i] = 0;
        bin_max[i] = 0;
      }
      bins_used = 0;
    }

    void print(){
      std:: cout << "Histogram:     "
                 << (EQ_W ? "eq_width" : "eq_height") << (block_wise ? "-blockwise" : "") << std:: endl
                 << "    Bins:      " << std::to_string(BINS_CNT) << std::endl
                 << "    Used bins: " << std::to_string(bins_used) << std::endl;

      auto j = 0;
      for(auto i = 0; i < BINS_CNT; i += j) {
        std::cout << "    ";
        for (j = 0; j < 8 && j+i < BINS_CNT; j++) {
          std::cout << "[" << std::to_string(bin_min[i + j]) << ";" << std::to_string(bin_max[i + j]) << ") ";
        }
        std::cout << std::endl;
      }

      std::cout << std::endl;
    }

    T get_min(){
      return bin_min[0];
    }

    T get_max(){
      return next_smaller_value(bin_max[bins_used-1]);
    }

    inline u64
    get_slot(const T v) const noexcept {
      //binary search
      $u32 low_bnd = 0;
      $u32 up_bnd = bins_used-1;
      $u32 slot = (up_bnd / 2);

      while(low_bnd != up_bnd){
        if(v < bin_max[slot]){
          if(v >= bin_min[slot])
            return slot;
          else{
            up_bnd--;
            slot = (up_bnd - low_bnd) / 2 + low_bnd;
          }
        }else{
          low_bnd++;
          slot = std::min((up_bnd-low_bnd)/2 + slot, bins_used-1);
        }
      }
      return slot;
    }

    u64 memory_footprint(){
      $u64 size_bins = BINS_CNT * sizeof(T) * 2;
      $u64 size_var = sizeof(bins_used);

      return size_bins + size_var;
    }
  };
};