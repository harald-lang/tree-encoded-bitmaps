#pragma once

#include <string>
#include <cstring>

class rule_of_five {

  char* cstring; // raw pointer used as a handle to a dynamically-allocated memory block

 public:

  explicit
  rule_of_five(const char* arg)
      : cstring(new char[std::strlen(arg) + 1]) {
    std::strcpy(cstring, arg); // populate
  }

  ~rule_of_five() {
    delete[] cstring;  // deallocate
  }

  rule_of_five(const rule_of_five& other) {
    cstring = new char[std::strlen(other.cstring) + 1];
    std::strcpy(cstring, other.cstring);
  }

  rule_of_five(rule_of_five&& other) : cstring(other.cstring) {
    other.cstring = nullptr;
  }

  rule_of_five&
  operator=(const rule_of_five& other) {
    char* tmp_cstring = new char[std::strlen(other.cstring) + 1];
    std::strcpy(tmp_cstring, other.cstring);
    delete[] cstring;
    cstring = tmp_cstring;
    return *this;
  }

  rule_of_five&
  operator=(rule_of_five&& other) {
    if (this != &other) {
      // prevent self-move
      delete[] cstring;
      cstring = other.cstring;
      other.cstring = nullptr;
    }
    return *this;
  }

};
