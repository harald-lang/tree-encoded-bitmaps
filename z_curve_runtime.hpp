#include "column_container.hpp"
using namespace std;

namespace z_curve_runtime{

//naive version
  template<typename T>
  bool cmp_coord_helper(const vector<cont::cont_elem>& tpl1,
                        const vector<cont::cont_elem>& tpl2,
                        const std::vector<$u32>& attr_idx_cmp) {
    $u16 j = 0;
    $u64 x = 0;

    $u16 numAttr = min(tpl1.size(), tpl2.size());

    for (auto i : attr_idx_cmp) { //determine the data type for padding
      T t1 = 0;
      T t2 = 0;
      switch(tpl1[i].size){
        case 8:
          t1 = *reinterpret_cast<$u8*>(tpl1[i].ptr);
          t2 = *reinterpret_cast<$u8*>(tpl2[i].ptr);
          break;
        case 16:
          t1 = *reinterpret_cast<$u16*>(tpl1[i].ptr);
          t2 = *reinterpret_cast<$u16*>(tpl2[i].ptr);
          break;
        case 32:
          t1 = *reinterpret_cast<$u32*>(tpl1[i].ptr);
          t2 = *reinterpret_cast<$u32*>(tpl2[i].ptr);
          break;
        case 64:
          t1 = *reinterpret_cast<$u64*>(tpl1[i].ptr);
          t2 = *reinterpret_cast<$u64*>(tpl2[i].ptr);
          break;
        default:
          cout << "Unsupported data type in cmp_coord_helper" << endl;
          exit(1);
      }

      auto y = t1 ^ t2;
      //cout << "Tuple 1: " << to_string(t1) << "      Tuple 2: " << to_string(t2) << "      XOR-Val: " << y << endl;

      if (x < y && (x < (x ^ y))) {
        j = i;
        x = y;
      }

    }

    $i64 t1 = 0;
    $i64 t2 = 0;

    switch(tpl1[j].size){
      case 8:
        t1 = *reinterpret_cast<$u8*>(tpl1[j].ptr);
        t2 = *reinterpret_cast<$u8*>(tpl2[j].ptr);
        break;
      case 16:
        t1 = *reinterpret_cast<$u16*>(tpl1[j].ptr);
        t2 = *reinterpret_cast<$u16*>(tpl2[j].ptr);
        break;
      case 32:
        t1 = *reinterpret_cast<$u32*>(tpl1[j].ptr);
        t2 = *reinterpret_cast<$u32*>(tpl2[j].ptr);
        break;
      case 64:
        t1 = *reinterpret_cast<$u64*>(tpl1[j].ptr);
        t2 = *reinterpret_cast<$u64*>(tpl2[j].ptr);
        break;
    }
    /*
    cout << "Tuple 1: " << to_string(t1) << "      Tuple 2: " << to_string(t2) << "      j: " << j <<
         "   Diff: " << (t1 - t2) <<"   Ret: "<< ((t1 - t2) < 0) << endl;
    */
    return ((t1 - t2) < 0);
  }

  bool calc_z_curve(const vector<cont::cont_elem>& tpl1,
                    const vector<cont::cont_elem>& tpl2,
                    const std::vector<$u32>& attr_idx_cmp) {
    $u8 max_col_size = 0;

    $u16 numAttr = min(tpl1.size(),tpl2.size());

//    for ($u16 i = 0; i < numAttr; i++) { //determine the data type for padding
    for (auto i : attr_idx_cmp) { //determine the data type for padding
      $u8 t1_type = tpl1[i].size;
      $u8 t2_type = tpl2[i].size;

      if (max(t1_type,t2_type) > max_col_size)
        max_col_size = t1_type;
    }

    bool cmp_result;

    switch (max_col_size) {
      case 8:
        cmp_result = cmp_coord_helper<$u8>(tpl1,tpl2, attr_idx_cmp);
        break;
      case 16:
        cmp_result = cmp_coord_helper<$u16>(tpl1,tpl2, attr_idx_cmp);
        break;
      case 32:
        cmp_result = cmp_coord_helper<$u32>(tpl1,tpl2, attr_idx_cmp);
        break;
      case 64:
        cmp_result = cmp_coord_helper<$u64>(tpl1,tpl2, attr_idx_cmp);
        break;
      default:
        cout << "Error, datatype not supported! (in calc_z_curve)" << endl;
        exit(1);
    }

    //cout << "Compare_result: " << cmp_result << endl;
    return cmp_result;
  }
// end naive version

  struct Compare{
    std::vector<$u32> attr_idx_cmp;
    bool operator()(const vector<cont::cont_elem>& tpl1, const vector<cont::cont_elem>& tpl2) const {
      return calc_z_curve(tpl1,tpl2, attr_idx_cmp);
    }
  };

  template <typename T>
  bool cmp_coord_helper(const vector<cont::container*>& tuples, u64 first_tuple, u64 second_tuple){

    if(first_tuple == second_tuple) {
      cout << "Just compare different tuples! (first_tuple != second_tuple)!" << endl;
      exit(1);
    }

    $u16 j = 0;
    $u64 x = 0;

    for($u16 i = 0; i < tuples.size(); i++) {

      T t1 = (T) (tuples[i]->get(first_tuple).data); //@Harald: (T) das sollte doch als padding reichen, oder?
      T t2 = (T) (tuples[i]->get(second_tuple).data);

      auto y = *t1 ^ *t2;
      //cout << "Tuple 1: " << to_string(*t1) << "      Tuple 2: " << to_string(*t2) << "      XOR-Val: " << y << endl;

      if(x<y && (x < (x ^ y))){
        j = i;
        x = y;
      }
    }

    T t1 = (T) (tuples[j]->get(first_tuple).data);
    T t2 = (T) (tuples[j]->get(second_tuple).data);

    return ((*t1 - *t2) < 0);
  }

  bool calc_z_curve(const vector<cont::container*>& tuples, u64 first_tuple, u64 second_tuple) {
    $u8 max_col_size = 0;

    for ($u8 i = 0; i < tuples.size(); i++) { //determine the data type for padding
      pair<bool, $u8> t1_type = tuples[i]->dataType();

      if (t1_type.second > max_col_size)
        max_col_size = t1_type.second;
    }

    bool cmp_result;

    switch (max_col_size) {
      case 8:
        cmp_result = cmp_coord_helper<$u8*>(tuples, first_tuple, second_tuple);
        break;
      case 16:
        cmp_result = cmp_coord_helper<$u16*>(tuples, first_tuple, second_tuple);
        break;
      case 32:
      cmp_result = cmp_coord_helper<$u32*>(tuples, first_tuple, second_tuple);
        break;
      case 64:
        cmp_result = cmp_coord_helper<$u64*>(tuples, first_tuple, second_tuple);
        break;
      default:
        cout << "Error, datatype not supported!" << endl;
        exit(1);
    }

    //cout << "Compare_result: " << cmp_result << endl;
    return cmp_result;
  }

  //old Version, which differs between signed and unsigned values (backup)
  bool calc_z_curve_diff_signed(const vector<cont::container*>& tuples, u64 first_tuple, u64 second_tuple) {

    bool is_signed = false;
    $u8 max_col_size = 0;

    for ($u8 i = 0; i < tuples.size(); i++) { //determine the data type for padding
      pair<bool, $u8> t1_type = tuples[i]->dataType();

      is_signed = is_signed || t1_type.first;

      if (t1_type.second > max_col_size)
        max_col_size = t1_type.second;
    }

    bool cmp_result;

    if (is_signed) {
      switch (max_col_size) {
        case 8:
          cmp_result = cmp_coord_helper<$i8*>(tuples, first_tuple, second_tuple);
          break;
        case 16:
          cmp_result = cmp_coord_helper<$i16*>(tuples, first_tuple, second_tuple);
          break;
        case 32:
          cmp_result = cmp_coord_helper<$i32*>(tuples, first_tuple, second_tuple);
          break;
        case 64:
          cmp_result = cmp_coord_helper<$i64*>(tuples, first_tuple, second_tuple);
          break;
        default:
          cout << "Error, datatype not supported!" << endl;
          exit(1);
      }
    } else {
      switch (max_col_size) {
        case 8:
          cmp_result = cmp_coord_helper<$u8*>(tuples, first_tuple, second_tuple);
          break;
        case 16:
          cmp_result = cmp_coord_helper<$u16*>(tuples, first_tuple, second_tuple);
          break;
        case 32:
          cmp_result = cmp_coord_helper<$u32*>(tuples, first_tuple, second_tuple);
          break;
        case 64:
          cmp_result = cmp_coord_helper<$u64*>(tuples, first_tuple, second_tuple);
          break;
        default:
          cout << "Error, datatype not supported!" << endl;
          exit(1);
      }
    }
    return cmp_result;
  }

};