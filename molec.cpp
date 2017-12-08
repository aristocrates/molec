#include "molec.hpp"
#include <iostream>

#define TEST_GAUSS 0

Point make_point(double x, double y) {
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

  //cout << box.get_nparticles() << endl;

  box.step();

  twodouble a;
  a.x = 1.;
  a.y = 1.;
  twodouble b;
  b.x = 3.;
  b.y = 3.;
  twodouble accel = box.acceleration(a, b);
  cout << "x: " << accel.x << ", y: " << accel.y << endl;

  twodouble c, d;
  c.x = 1.;
  c.y = 10.;
  d.x = 29.;
  d.y = 10.;
  accel = box.acceleration(c, d);
  cout << "x: " << accel.x << ", y: " << accel.y << endl;
  accel = box.acceleration(d, c);
  cout << "x: " << accel.x << ", y: " << accel.y << endl;

  c.x = 5.;
  c.y = 29.;
  d.x = 5.;
  d.y = 1.;
  accel = box.acceleration(c, d);
  cout << "x: " << accel.x << ", y: " << accel.y << endl;
  accel = box.acceleration(d, c);
  cout << "x: " << accel.x << ", y: " << accel.y << endl;

  cout << "Two point system test" << endl;

  Point centers2[2] = {make_point(10., 10.), make_point(10., 11.)};
  LennardJonesSystem box2(2, centers2, 1.0, 2.0, 0.1, 0., 0.2, 0.001, 30);
  cout << "(" << box2.get_centers()[0].x << "," << box2.get_centers()[0].y << ")";
  cout << "(" << box2.get_centers()[1].x << "," << box2.get_centers()[1].y << ")" << endl;

  for (int i = 0; i < 200; i++) {
    box2.step();
  }

  cout << "(" << box2.get_centers()[0].x << "," << box2.get_centers()[0].y << ")";
								    cout << "(" << box2.get_centers()[1].x << "," << box2.get_centers()[1].y << ")" << endl;

  cout << "Two point system sanity test" << endl;

  double centers3[4] = {10., 10., 10., 11.};
  LennardJonesSystem box3(2, centers2, 1.0, 2.0, 0.1, 0., 0.2, 0.001, 30);
  cout << "(" << box3.get_centers()[0].x << "," << box2.get_centers()[0].y << ")";
  cout << "(" << box3.get_centers()[1].x << "," << box2.get_centers()[1].y << ")" << endl;

  for (int i = 0; i < 200; i++) {
    box3.step();
  }

  cout << "(" << box3.get_centers()[0].x << "," << box3.get_centers()[0].y << ")";
  cout << "(" << box3.get_centers()[1].x << "," << box3.get_centers()[1].y << ")" << endl;

  cout << "Testing three point system" << endl;



  #if TEST_GAUSS
  cout << "Testing Gaussian random noise" << endl;

  cout << "[";
  for (int i = 0; i < 1000; i++) {
    twodouble gaussian_rand = normal_gaussian_rand();
    cout << gaussian_rand.x << ", " << gaussian_rand.y << ", ";
  }
  twodouble gaussian_rand = normal_gaussian_rand();
  cout << gaussian_rand.x << ", " << gaussian_rand.y << "]" << endl;
  #endif
  return 0;
}
