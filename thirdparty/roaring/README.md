# CRoaringUnityBuild
Dumps of CRoaring unity builds (for convenience)

This code is automatically generated from https://github.com/RoaringBitmap/CRoaring


## Usage (C)

```C
#include <stdio.h>
#include "roaring.c"
int main() {
  roaring_bitmap_t *r1 = roaring_bitmap_create();
  for (uint32_t i = 100; i < 1000; i++) roaring_bitmap_add(r1, i);
  printf("cardinality = %d\n", (int) roaring_bitmap_get_cardinality(r1));
  roaring_bitmap_free(r1);
  return 0;
}
```
## Usage (C++)


```C++
#include <iostream>
#include "roaring.hh"
#include "roaring.c"
int main() {
  Roaring r1;
  for (uint32_t i = 100; i < 1000; i++) {
    r1.add(i);
  }
  std::cout << "cardinality = " << r1.cardinality() << std::endl;
  return 0;
}

