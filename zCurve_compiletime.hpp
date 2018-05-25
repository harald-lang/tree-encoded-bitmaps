#include "dtl/adept.hpp"
#include <assert.h>
#include <tuple>

using namespace std;

//abgewandelt von: https://stackoverflow.com/questions/28997271/c11-way-to-index-tuple-at-runtime-without-using-switch
template<size_t I>
struct visit_impl {
  template<typename T, typename F>
  static $u16 visit(T& tup, size_t idx, F fun) {
    if (idx == I - 1) fun(std::get<I - 1>(tup));
    else visit_impl<I - 1>::visit(tup, idx, fun);
  }
};

template<>
struct visit_impl<0> {
  template<typename T, typename F>
  static $u16 visit(T& tup, size_t idx, F fun) { assert(false); }
};

template<typename F, typename... Ts>
$u16 visit_at(std::tuple<Ts...> const& tup, size_t idx, F fun) {
  visit_impl<sizeof...(Ts)>::visit(tup, idx, fun);
}

template<typename F, typename... Ts>
$u16 visit_at(std::tuple<Ts...>& tup, size_t idx, F fun) {
  visit_impl<sizeof...(Ts)>::visit(tup, idx, fun);
}

template <typename... T>
bool cmp_z_coord(const tuple<T...>& tup1, const tuple<T...>& tup2){

  auto attr_number = sizeof...(T);
  $u16 j = 0;
  $u64 x = 0;

  auto get_elem = [](auto t1) {return t1;};

  for($u16 i = 0; i < attr_number; i++){
    auto y = visit_at(tup1,i,get_elem) ^ visit_at(tup2,i,get_elem);

    //cout << "tup1: " << visit_at(tup1,i,get_elem) << " | tup2: " << visit_at(tup2,i,get_elem) << " | y: " << y;

    if(x<y && (x < (x ^ y))){
      //cout << " | New j:" << i ;
      j = i;
      x = y;
    }
    //cout << endl;
  }
  //cout << "Final j: " << j << endl;
  //cout << "Difference: " << (visit_at(tup1,j,get_elem) - visit_at(tup2,j,get_elem)) << endl;

  return ((visit_at(tup1,j,get_elem) - visit_at(tup2,j,get_elem)) < 0);
}

struct Compare{
  template<typename... T>
  bool operator()(const tuple<T...>& tup1, const tuple<T...>& tup2) const { return cmp_z_coord(tup1,tup2);}
};



