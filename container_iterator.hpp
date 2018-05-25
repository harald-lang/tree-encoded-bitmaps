#pragma once

#include <dtl/adept.hpp>
#include <iterator>
#include "column_container.hpp"


class cont_iter : public std::iterator<random_access_iterator_tag, cont::cont_elem>{

  vector<cont::cont_elem> columns;

  public:
    using difference_type = typename std::iterator<std::random_access_iterator_tag, cont::cont_elem>::difference_type;

    cont_iter(vector<cont::container> src){
      for($u16 i = 0; i < src.size(); i++) {
        columns.push_back(cont::cont_elem(src[i].getParam()));
      }
    }

    cont_iter(const cont_iter &rhs){
      for($u16 i = 0; i < rhs.columns.size(); i++) {
        columns.push_back(cont::cont_elem(rhs.columns[i].size, rhs.columns[i].ptr));
      }
    }

    /* inline Iterator& operator=(Type* rhs) {_ptr = rhs; return *this;} */
    /* inline Iterator& operator=(const Iterator &rhs) {_ptr = rhs._ptr; return *this;} */

    cont_iter& operator+=(difference_type rhs) {

      for(int i = 0; i < columns.size(); i++){
        switch(columns[i].size){
          case 8:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + (1*rhs);
            break;
          case 16:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + (2*rhs);
            break;
          case 32:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + (4*rhs);
            break;
          case 64:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + (8*rhs);
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }

      return *this;
    }

    cont_iter& operator-=(difference_type rhs) {

      for(int i = 0; i < columns.size(); i++){
        switch(columns[i].size){
          case 8:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - (1*rhs);
            break;
          case 16:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - (2*rhs);
            break;
          case 32:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - (4*rhs);
            break;
          case 64:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - (8*rhs);
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }

      return *this;
    }


    //inline cont::cont_elem& operator*() const {return *_ptr;}
    //inline cont::cont_elem* operator->() const {return _ptr;}
    vector<cont::cont_elem> operator[](difference_type rhs) const {
      vector<cont::cont_elem> tmp;

      for(int i = 0; i < columns.size(); i++){
        switch(columns[i].size){
          case 8:
            tmp.push_back(cont::cont_elem(8,static_cast<$u8*>(columns[i].ptr) + rhs));
            break;
          case 16:
            tmp.push_back(cont::cont_elem(16,static_cast<$u16*>(columns[i].ptr) + rhs));
            break;
          case 32:
            tmp.push_back(cont::cont_elem(32,static_cast<$u32*>(columns[i].ptr) + rhs));
            break;
          case 64:
            tmp.push_back(cont::cont_elem(64,static_cast<$u64*>(columns[i].ptr) + rhs));
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }

      return tmp;
      //return _ptr[rhs];
    }

    cont_iter& operator++() {

      for(int i = 0; i < columns.size(); i++){
        switch(columns[i].size){
          case 8:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + 1;
            break;
          case 16:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + 2;
            break;
          case 32:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + 4;
            break;
          case 64:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) + 8;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return *this;
    }

    cont_iter& operator--() {

      for(int i = 0; i < columns.size(); i++){
        switch(columns[i].size){
          case 8:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - 1;
            break;
          case 16:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - 2;
            break;
          case 32:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - 4;
            break;
          case 64:
            columns[i].ptr = static_cast<char*>(columns[i].ptr) - 8;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }

      return *this;
    }

    //cont_iter operator++(int) const {cont_iter tmp(*this); ++_ptr; return tmp;}

    //inline cont_iter operator--(int) const {cont_iter tmp(*this); --_ptr; return tmp;} //?

    //inline cont_iter operator+(const cont_iter& rhs) {return cont_iter(_ptr+rhs.ptr);}

    difference_type operator-(const cont_iter& rhs) const {
      long n = 0;

      switch(columns[0].size){
        case 8:
          n = static_cast<$u8*>(rhs.columns[0].ptr) - static_cast<$u8*>(columns[0].ptr);
          break;
        case 16:
          n = static_cast<$u16*>(rhs.columns[0].ptr) - static_cast<$u16*>(columns[0].ptr);
          break;
        case 32:
          n = static_cast<$u32*>(rhs.columns[0].ptr) - static_cast<$u32*>(columns[0].ptr);
          break;
        case 64:
          n = static_cast<$u64*>(rhs.columns[0].ptr) - static_cast<$u64*>(columns[0].ptr);
          break;
        default:
          cout << "Error, datatype not supported!" << endl;
          exit(1);
      }

      return n;
    }

    //inline cont_iter operator+(difference_type rhs) const { return cont_iter(_ptr+rhs); }
    //inline cont_iter operator-(difference_type rhs) const {return cont_iter(_ptr-rhs);}
    //friend inline cont_iter operator+(difference_type lhs, const cont_iter& rhs) {return cont_iter(lhs+rhs._ptr);}
    //friend inline cont_iter operator-(difference_type lhs, const cont_iter& rhs) {return cont_iter(lhs-rhs._ptr);}

    bool operator==(const cont_iter& rhs) const {

      if(columns.size() != rhs.columns.size()) //check #columns
        return false;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return false;

        switch(columns[i].size) { //check all columns
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) != static_cast<$u8*>(columns[0].ptr))
              return false;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) != static_cast<$u16*>(columns[0].ptr))
              return false;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) != static_cast<$u32*>(columns[0].ptr))
              return false;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) != static_cast<$u64*>(columns[0].ptr))
              return false;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return true;
    }

    bool operator!=(const cont_iter& rhs) const {

      if(columns.size() != rhs.columns.size()) //check #columns
        return true;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return true;

        switch(columns[i].size) {
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) != static_cast<$u8*>(columns[0].ptr))
              return true;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) != static_cast<$u16*>(columns[0].ptr))
              return true;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) != static_cast<$u32*>(columns[0].ptr))
              return true;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) != static_cast<$u64*>(columns[0].ptr))
              return true;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return false;
    }

    bool operator>(const cont_iter& rhs) const {
      if(columns.size() != rhs.columns.size()) //check #columns
        return false;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return false;

        switch(columns[i].size) { //check all columns
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) >= static_cast<$u8*>(columns[0].ptr))
              return false;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) >= static_cast<$u16*>(columns[0].ptr))
              return false;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) >= static_cast<$u32*>(columns[0].ptr))
              return false;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) >= static_cast<$u64*>(columns[0].ptr))
              return false;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return true;
    }

    bool operator<(const cont_iter& rhs) const {
      if(columns.size() != rhs.columns.size()) //check #columns
        return false;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return false;

        switch(columns[i].size) { //check all columns
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) <= static_cast<$u8*>(columns[0].ptr))
              return false;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) <= static_cast<$u16*>(columns[0].ptr))
              return false;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) <= static_cast<$u32*>(columns[0].ptr))
              return false;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) <= static_cast<$u64*>(columns[0].ptr))
              return false;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return true;
    }

    bool operator>=(const cont_iter& rhs) const { //_ptr >= rhs._ptr;

      if(columns.size() != rhs.columns.size()) //check #columns
        return false;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return false;

        switch(columns[i].size) { //check all columns
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) > static_cast<$u8*>(columns[0].ptr))
              return false;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) > static_cast<$u16*>(columns[0].ptr))
              return false;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) > static_cast<$u32*>(columns[0].ptr))
              return false;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) > static_cast<$u64*>(columns[0].ptr))
              return false;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return true;}

    inline bool operator<=(const cont_iter& rhs) const { //_ptr <= rhs._ptr;
      if(columns.size() != rhs.columns.size()) //check #columns
        return false;

      for(int i = 0; i < columns.size(); i++){
        if(columns[i].size != rhs.columns[i].size) //check data type
          return false;

        switch(columns[i].size) { //check all columns
          case 8:
            if(static_cast<$u8*>(rhs.columns[0].ptr) < static_cast<$u8*>(columns[0].ptr))
              return false;
            break;
          case 16:
            if(static_cast<$u16*>(rhs.columns[0].ptr) < static_cast<$u16*>(columns[0].ptr))
              return false;
            break;
          case 32:
            if(static_cast<$u32*>(rhs.columns[0].ptr) < static_cast<$u32*>(columns[0].ptr))
              return false;
            break;
          case 64:
            if(static_cast<$u64*>(rhs.columns[0].ptr) < static_cast<$u64*>(columns[0].ptr))
              return false;
            break;
          default:
            cout << "Error, datatype not supported!" << endl;
            exit(1);
        }
      }
      return true;
    }
};