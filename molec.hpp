/**
 * 
 */
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <cassert>

/**
 * Encapsulates a Langevin integrator
 */
class Langevin {
public:
  Langevin(float T, float alpha);
  ~Langevin();
private:
  float T, alpha;
};

Langevin::Langevin(float T, float alpha) {
  assert(T > 0);
  this->T = T;
  this->alpha = alpha;
}

// The first dimension will be y, the second will be x
#define get(x, y, n, grid) grid[y % n][x % n];
/**
 * Holds the value of the field at a series of discrete points (in an nxn grid)
 * with periodic boundary conditions.
 */
class Mesh {
public:
  Mesh(int n, double delta);
  ~Mesh();
private:
  int n;
  float delta;
  float **grid;
};

Mesh::Mesh(int n, double delta) {
  assert(n > 0);
  this->n = n;
  this->delta = delta;
  grid = new float*[n];
  for (int i = 0; i < n; i++) {
    grid[i] = new float[n];
  }
}

Mesh::~Mesh() {
  for (int i = 0; i < n; i++) {
    delete[] grid[i];
  }
  delete[] grid;
}

float 

#endif
