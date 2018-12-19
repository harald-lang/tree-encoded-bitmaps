#pragma once

#include <dtl/adept.hpp>
#include <memory>

//Covered datatypes: $u8, $u16, $u32, $u64, $i8, $i16, $i32, $i64
namespace cont{

struct container_format{
  char* data;
  bool is_signed;
  $u8 size_bits;

  container_format(char* data_ptr, bool is_Signed, $u8 size_in_bits){
    data = data_ptr;
    is_signed = is_Signed;
    size_bits = size_in_bits;
  }
};

struct cont_elem{
  $u8 size;
  char* ptr;

  cont_elem($u8 bits, char* tpl){
    size = bits;
    ptr = tpl;
  }

  void print(){
    switch(size){
      case 8:
        std::cout << std::to_string(*(reinterpret_cast<$u8*>(ptr)));
        break;
      case 16:
        std::cout << std::to_string(*(reinterpret_cast<$u16*>(ptr)));
        break;
      case 32:
        std::cout << std::to_string(*(reinterpret_cast<$u32*>(ptr)));
        break;
      case 64:
        std::cout << std::to_string(*(reinterpret_cast<$u64*>(ptr)));
        break;
      default:
        std::cout << "Unsupported data type! (in cont_elem.print())";
    }
  }
};

struct container{
  container(){};
  virtual void print_Data()=0;
  virtual container_format get(u32 index) = 0;
  virtual std::pair<bool, u8> dataType()=0;
  virtual cont_elem getParam() = 0;
};

struct $u8_column : virtual container{
  $u8 *column;
  $u64 column_size;
public:
  $u8_column(void *col, $u64 num_tup) {
    column = ($u8*) col;
    column_size = num_tup;
  };

  std::pair<bool,u8> dataType() override {
    return {false,8}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), false, sizeof($u8));
  }

  cont_elem getParam() override {
    return cont_elem(8,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $u16_column : virtual container{
  $u16 *column;
  $u64 column_size;
public:
  $u16_column(void *col, $u64 num_tup) {
    column = ($u16*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {false,16}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), false, sizeof($u16));
  }


  cont_elem getParam() override {
    return cont_elem(16,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $u32_column : virtual container{
  $u32 *column;
  $u64 column_size;
public:
  $u32_column(void *col, $u64 num_tup) {
    column = ($u32*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {false,32}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), false, sizeof($u32));
  }

  cont_elem getParam() override {
    return cont_elem(32,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $u64_column : virtual container{
  $u64 *column;
  $u64 column_size;
public:
  $u64_column(void *col, $u64 num_tup) {
    column = ($u64*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {false, 64}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), false, sizeof($u64));
  }

  cont_elem getParam() override {
    return cont_elem(64,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $i8_column : virtual container{
  $i8 *column;
  $u64 column_size;
public:
  $i8_column(void *col, $u64 num_tup) {
    column = ($i8*) col;
    column_size = num_tup;
  };

  std::pair<bool,u8> dataType() override {
    return {true, 8}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), true, sizeof($i8));
  }

  cont_elem getParam() override {
    return cont_elem(8,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $i16_column : virtual container{
  $i16 *column;
  $u64 column_size;
public:
  $i16_column(void *col, $u64 num_tup) {
    column = ($i16*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {true, 16}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), true, sizeof($i16));
  }

  cont_elem getParam() override {
    return cont_elem(16,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $i32_column : virtual container{
  $i32 *column;
  $u64 column_size;
public:
  $i32_column(void *col, $u64 num_tup) {
    column = ($i32*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {true, 32}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), true, sizeof($i32));
  }


  cont_elem getParam() override {
    return cont_elem(32,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

struct $i64_column : virtual container{
  $i64 *column;
  $u64 column_size;
public:
  $i64_column(void *col, $u64 num_tup) {
    column = ($i64*) col;
    column_size = num_tup;
  }

  std::pair<bool,u8> dataType() override {
    return {true, 64}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), true, sizeof($i64));
  }


  cont_elem getParam() override {
    return cont_elem(64,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      std::cout << std::to_string(column[i]) << std::endl;
  }
};

/*
struct $u1_column : virtual container{
  $u1 *column;
  $u64 column_size;

  $u1_column(void *col, $u64 num_tup) {
    column = ($u1*) col;
    column_size = num_tup;
  }

  pair<bool, u8> dataType() override {
    return {false,1}; //returns
  };

  container_format get(u32 index) override {
    return container_format(reinterpret_cast<char*>(&column[index]), false, sizeof($u1));
  }

  cont_elem getParam() override {
    return cont_elem(1,reinterpret_cast<char*>(column));
  }

  void print_Data() override {
    for(int i = 0; i < column_size; i++)
      cout << to_string(column[i]) << endl;
  }
};
unique_ptr<container> storeColumn($u1 *col, $u64 num_elements){
  container* column = new $u1_column(col, num_elements);
  return unique_ptr<container>(column);
}
*/

std::unique_ptr<container> storeColumn($u8 *col, $u64 num_elements){
  container* column = new $u8_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($u16 *col, $u64 num_elements){
  container* column = new $u16_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($u32 *col, $u64 num_elements){
  container* column = new $u32_column(col, num_elements);

  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($u64 *col, $u64 num_elements){
  container* column = new $u64_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($i8 *col, $u64 num_elements){
  container* column = new $i8_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($i16 *col, $u64 num_elements){
  container* column = new $i16_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($i32 *col, $u64 num_elements){
  container* column = new $i32_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

std::unique_ptr<container> storeColumn($i64 *col, $u64 num_elements){
  container* column = new $i64_column(col, num_elements);
  return std::unique_ptr<container>(column);
}

//applies the permutation to an array for further tasks
void applyPermutation(std::vector<std::vector<cont::cont_elem>> rel_perm, std::vector<char*> new_rel){

  for(auto row = 0; row < rel_perm.size(); row++){
    for(auto att = 0; att < rel_perm[row].size(); att++){
      switch(rel_perm[row][att].size){
        case 8:
          new_rel[att][row] = *rel_perm[row][att].ptr;
          break;
        case 16:
          (reinterpret_cast<$u16*>(new_rel[att]))[row] = *reinterpret_cast<$u16*>(rel_perm[row][att].ptr);
          break;
        case 32:
          (reinterpret_cast<$u32*>(new_rel[att]))[row] = *reinterpret_cast<$u32*>(rel_perm[row][att].ptr);
          break;
        case 64:
          (reinterpret_cast<$u64*>(new_rel[att]))[row] = *reinterpret_cast<$u64*>(rel_perm[row][att].ptr);
          break;
        default:
          std::cout << "Unsupported data type! applyPermutation(...)" << std::endl;
      }
    }
  }
}

} // namespace cont