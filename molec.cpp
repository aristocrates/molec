#include "molec.hpp"
#include <iostream>

Point make_point(float x, float y) {
  Point ans;
  ans.x = x;
  ans.y = y;
  return ans;
}

int main() {
  Point centers[10] = {make_point(1, 2), make_point(1.4, 2.7),
		       make_point(4, 4), make_point(5, 5),
		       make_point(8, 8), make_point(9, 9),
		       make_point(10, 10), make_point(1, 5),
		       make_point(15, 10), make_point(10, 15)};
  LennardJonesSystem box(10, centers, 1.0, 1.0, 1.0, 100, 0.2, 0.001, 30);
  return 0;
}
