/**
 * 
 */
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <cassert>
#include <cmath>
#include <cstdlib>

// Constants
#define KBOLTZMANN (8.6 * 1e-5)

typedef struct {
  float x, y;
} twofloat;

typedef twofloat Point;

/**
 * Encapsulates a Langevin integrator for the Lennard Jones potential in
 * two dimensions.
 *
 * One of these is created for every particle.
 */
class Langevin {
public:
  Langevin();
  Langevin(float T, float alpha, float h, Point *p);
  ~Langevin();
  void param_change();
  twofloat update(twofloat acceleration);
private:
  Point *p;
  bool first;
  twofloat position, oldposition;
  twofloat prevrand;
  float T, alpha, h, time;
  float A, B, C, D;
};

/**
 * Default constructor just to allow for array creation
 */
Langevin::Langevin() {
  // do nothing at all
}

Langevin::Langevin(float T, float alpha, float h, Point *p) {
  assert(T > 0);
  this->T = T;
  this->alpha = alpha;
  this->h = h;

  param_change();

  this->time = 0.;

  this->prevrand.x = 0;
  this->prevrand.y = 0;
  this->first = true;
  this->p = p;
  this->oldposition.x = p->x;
  this->oldposition.y = p->y;
  this->position.x = p->x;
  this->position.y = p->y;
}

Langevin::~Langevin() {
}

/**
 * Whenever alpha, the timestep h, or T change, this must be run
 */
void Langevin::param_change() {
  this->A = 4.0 / (2.0 + this->alpha * this->h);
  this->B = -1.0 * (2.0 - this->alpha * this->h) / (2.0 + this->alpha * this->h);
  this->C = 2.0 * this->h * this->h / (2.0 + this->alpha * this->h);
  this->D = sqrt(2 * this->alpha * KBOLTZMANN * this->T) / (2.0 + this->alpha * this->h);
}

float random_uniform(float a, float b) {
  float zero_to_one = 1.0 * (random() % (1<<30)) / (1<<30);
  return a + zero_to_one * b;
}

/**
 * Applies the Box-Muller transformation to two uniform random numbers
 * to get two gaussian random numbers
 */
twofloat normal_gaussian_rand() {
  twofloat ans;

  float x1 = random_uniform(0, 1);
  float x2 = random_uniform(0, 1);

  ans.x = sqrt( -2.0 * log(x1) ) * cos(2 * M_PI * x2);
  ans.y = sqrt( -2.0 * log(x1) ) * sin(2 * M_PI * x2);

  return ans;
}

/**
 * Updates the point given acceleration
 */
twofloat Langevin::update(twofloat acceleration) {
  twofloat ans;

  twofloat current_rand = normal_gaussian_rand();

  float newx = (this->A * this->p->x + this->B * this->oldposition.x
		+ this->C * acceleration.x
		+ this->D * (current_rand.x + this->prevrand.x) );
  float newy = (this->A * this->p->y + this->B * this->oldposition.y
		+ this->C * acceleration.y
		+ this->D * (current_rand.y + this->prevrand.y) );

  this->oldposition.x = this->p->x;
  this->oldposition.y = this->p->y;

  this->p->x = newx;
  this->p->y = newy;

  prevrand = current_rand;
  return ans;
}

class LennardJonesSystem {
public:
  LennardJonesSystem(int nparticles, Point *centers, float m, float rmin,
		     float epsilon, float T, float alpha, float h,
		     float rcmult, float box_width);
  LennardJonesSystem(int nparticles, Point *centers, float m, float rmin,
		     float epsilon, float T, float alpha, float h,
		     float box_width);
  ~LennardJonesSystem();
  twofloat acceleration(twofloat from, twofloat on);
  void step();
private:
  float box_width;
  float T, alpha, h;
  float m, rmin, rc, epsilon;

  int nparticles;
  Point *centers;
  twofloat *accels;
  Langevin *integrators;
};

#define DEFAULT_RC_MULT 2.3

LennardJonesSystem::LennardJonesSystem(int nparticles, Point *centers,
				       float m, float rmin,
				       float epsilon, float T, float alpha,
				       float h, float box_width) {
  LennardJonesSystem(nparticles, centers, m, rmin, epsilon, T, alpha, h,
		     DEFAULT_RC_MULT, box_width);
}

LennardJonesSystem::LennardJonesSystem(int nparticles, Point *centers, float m,
				       float rmin, float epsilon,
				       float T, float alpha, float h,
				       float rcmult, float box_width) {
  this->box_width = box_width;

  this->m = m;
  this->rmin = rmin;
  this->rc = rmin * rcmult;
  this->epsilon = epsilon;

  this->T = T;
  this->alpha = alpha;
  this->h = h;

  this->nparticles = nparticles;
  this->centers = centers;
  this->accels = new twofloat[nparticles];
  for (int i = 0; i < nparticles; i++) {
    this->accels[i].x = 0;
    this->accels[i].y = 0;
  }
  this->integrators = new Langevin[nparticles];
  for (int i = 0; i < nparticles; i++) {
    this->integrators[i] = Langevin(this->T, this->alpha, this->h,
				    &centers[i]);
  }
}

LennardJonesSystem::~LennardJonesSystem() {
  delete[] this->accels;
  delete[] this->integrators;
}

/**
 * Calculates the component of acceleration at point "on" due to the
 * contribution of the "from" point.
 *
 * The following antisymmetry holds: acceleration(a, b) = -acceleration(b, a)
 */
twofloat LennardJonesSystem::acceleration(Point from, Point on) {
  float dx = on.x - from.x;
  float dy = on.y - from.y;
  if (dx > this->box_width / 2) {
    dx -= this->box_width;
  }
  if (dy > this->box_width / 2) {
    dy -= this->box_width;
  }
  float r = sqrt(dx * dx + dy * dy);
  twofloat ans;
  if (r > rc) {
    ans.x = 0;
    ans.y = 0;
    return ans;
  }
  float magnitude = pow(rmin, 6) / pow(r, 6) - 1.0;
  magnitude *= 12.0 * this->epsilon / this->m * pow(rmin, 6) / pow(r, 7);
  ans.x = magnitude * dx / r;
  ans.y = magnitude * dy / r;
  return ans;
}

#endif
