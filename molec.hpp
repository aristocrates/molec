/**
 * Molecular dynamics simulator
 */
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

// https://stackoverflow.com/questions/11980292/how-to-wrap-around-a-range
#define REALMOD(x, y) ( (x) - (y) * floor( (x) / (y) ) )

// Constants
#define KBOLTZMANN (1.0)

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
  Langevin(float T, float alpha, float h, float box_width, Point *p);
  ~Langevin();
  void param_change(float T, float alpha, float h, float box_width);
  void update(twofloat acceleration);
private:
  Point *p;
  //bool first;
  twofloat position, oldposition;
  twofloat prevrand;
  float T, alpha, h, box_width, time;
  float A, B, C, D;
};

/**
 * Default constructor just to allow for array creation
 */
Langevin::Langevin() {
  // do nothing at all
}

Langevin::Langevin(float T, float alpha, float h, float box_width, Point *p) {
  assert(T >= 0);
  assert(alpha > 0);
  assert(h > 0);

  param_change(T, alpha, h, box_width);

  this->time = 0.;

  this->prevrand.x = 0;
  this->prevrand.y = 0;
  //this->first = true;
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
void Langevin::param_change(float T, float alpha, float h, float box_width) {
  this->T = T;
  this->alpha = alpha;
  this->h = h;
  this->box_width = box_width;

  this->A = 4.0 / (2.0 + alpha * h);
  this->B = -1.0 * (2.0 - alpha * h) / (2.0 + alpha * h);
  this->C = 2.0 * h * h / (2.0 + alpha * h);
  this->D = sqrt(2 * alpha * KBOLTZMANN * T * h * h * h) / (2.0 + alpha * h);
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
void Langevin::update(twofloat acceleration) {
  twofloat current_rand = normal_gaussian_rand();

  float newx = (this->A * this->p->x + this->B * this->oldposition.x
		+ this->C * acceleration.x
		+ this->D * (current_rand.x + this->prevrand.x) );
  float newy = (this->A * this->p->y + this->B * this->oldposition.y
		+ this->C * acceleration.y
		+ this->D * (current_rand.y + this->prevrand.y) );

  // If newx or newy need to wrap around, then correct for that
  float newx_mod = REALMOD(newx, this->box_width);
  float newy_mod = REALMOD(newy, this->box_width);
  // Also need to correct the previous location
  // Failing to do this will give the particle a sudden "kick" acceleration
  float jumpx = newx_mod - newx;
  float jumpy = newy_mod - newy;
  // manually set the jumps to zero if they are just noise
  if (abs(jumpx) < this->box_width / 2) {
    jumpx = 0;
  }
  if (abs(jumpy) < this->box_width / 2) {
    jumpy = 0;
  }

  this->oldposition.x = this->p->x + jumpx;
  this->oldposition.y = this->p->y + jumpy;

  this->p->x = newx_mod;
  this->p->y = newy_mod;

  this->prevrand = current_rand;
}

class LennardJonesSystem {
public:
  LennardJonesSystem(int nparticles, Point *centers, float m, float rmin,
		     float epsilon, float T, float alpha, float h,
		     float rc, float box_width);
  LennardJonesSystem(int nparticles, Point *centers, float m, float rmin,
		     float epsilon, float T, float alpha, float h,
		     float box_width);
  LennardJonesSystem(int nparticles, float *centers, float m,
		     float rmin, float epsilon, float T,
		     float alpha, float h, float rc,
		     float box_width);
  LennardJonesSystem(int nparticles, float *centers, float m,
		     float rmin, float epsilon, float T,
		     float alpha, float h, float box_width);
  ~LennardJonesSystem();
  twofloat acceleration(twofloat from, twofloat on);
  void param_change(float T, float alpha, float h, float box_width,
		    float m, float rmin, float rc, float epsilon);
  void param_change(float T, float alpha, float h, float box_width);
  void step();
  float get_box_width() { return this->box_width;}
  float get_T() { return this->T;}
  float get_alpha() { return this->alpha;}
  float get_h() { return this->h;}
  float get_m() { return this->m;}
  float get_rmin() { return this->rmin;}
  float get_rc() { return this->rc;}
  float get_epsilon() { return this->epsilon;}
  float get_nparticles() { return this->nparticles;}
  twofloat *get_centers() { return this->centers;}
private:
  float box_width;
  float T, alpha, h;
  float m, rmin, rc, epsilon;
  bool delete_centers = false;

  int nparticles;
  Point *centers;
  twofloat *accels;
  Langevin *integrators;

  void construct_lennard_jones(int nparticles, Point *centers, float m,
			       float rmin, float epsilon, float T, float alpha,
			       float h, float rc, float box_width);
  void param_set(float T, float alpha, float h, float box_width,
		 float m, float rmin, float rc, float epsilon);
  void update_integrators();
};

#define DEFAULT_RC_MULT 2.3

LennardJonesSystem::LennardJonesSystem(int nparticles, Point *centers,
				       float m, float rmin,
				       float epsilon, float T, float alpha,
				       float h, float box_width) {
  this->delete_centers = false;
  construct_lennard_jones(nparticles, centers, m, rmin, epsilon,
			  T, alpha, h, DEFAULT_RC_MULT * rmin, box_width);
}

LennardJonesSystem::LennardJonesSystem(int nparticles, Point *centers, float m,
				       float rmin, float epsilon,
				       float T, float alpha, float h,
				       float rc, float box_width) {
  this->delete_centers = false;
  construct_lennard_jones(nparticles, centers, m, rmin, epsilon,
			  T, alpha, h, rc, box_width);
}

/**
 * Constructors accepting an array of floats for centers, makes using ctypes
 * easier.
 */
LennardJonesSystem::LennardJonesSystem(int nparticles, float *centers, float m,
				       float rmin, float epsilon, float T,
				       float alpha, float h, float rc,
				       float box_width) {
  this->delete_centers = true;
  Point *local_centers = new Point[nparticles];
  for (int i = 0; i < nparticles; i++) {
    local_centers[i].x = centers[2 * i];
    local_centers[i].y = centers[2 * i + 1];
  }
  construct_lennard_jones(nparticles, local_centers, m, rmin, epsilon, T,
			  alpha, h, rc, box_width);
}
LennardJonesSystem::LennardJonesSystem(int nparticles, float *centers, float m,
				       float rmin, float epsilon, float T,
				       float alpha, float h, float box_width) {
  this->delete_centers = true;
  Point *local_centers = new Point[nparticles];
  for (int i = 0; i < nparticles; i++) {
    local_centers[i].x = centers[2 * i];
    local_centers[i].y = centers[2 * i + 1];
  }
  construct_lennard_jones(nparticles, local_centers, m, rmin, epsilon, T,
			  alpha, h, DEFAULT_RC_MULT * rmin, box_width);
}

void LennardJonesSystem::construct_lennard_jones(int nparticles,
						 Point *centers, float m,
						 float rmin, float epsilon,
						 float T, float alpha, float h,
						 float rc, float box_width) {
  param_set(T, alpha, h, box_width, m, rmin, rc, epsilon);

  this->nparticles = nparticles;
  //cout << this->nparticles << endl;
  this->centers = centers;
  this->accels = new twofloat[nparticles];
  for (int i = 0; i < nparticles; i++) {
    this->accels[i].x = 0;
    this->accels[i].y = 0;
  }
  //cout << this->nparticles << endl;
  this->integrators = new Langevin[nparticles];
  for (int i = 0; i < nparticles; i++) {
    this->integrators[i] = Langevin(this->T, this->alpha, this->h,
				    this->box_width, &centers[i]);
  }
  //cout << this->nparticles << endl;
}

LennardJonesSystem::~LennardJonesSystem() {
  if (delete_centers) {
    delete[] this->centers;
  }
  delete[] this->accels;
  delete[] this->integrators;
}

void LennardJonesSystem::param_change(float T, float alpha, float h,
				      float box_width) {
  this->T = T;
  this->alpha = alpha;
  this->h = h;
  this->box_width = box_width;

  // Update all the integrators
  update_integrators();
}


void LennardJonesSystem::param_change(float T, float alpha, float h,
				      float box_width, float m, float rmin,
				      float rc, float epsilon) {
  param_set(T, alpha, h, box_width, m, rmin, rc, epsilon);

  // Update all the integrators
  update_integrators();
}

void LennardJonesSystem::update_integrators() {
  for (int i = 0; i < this->nparticles; i++) {
    this->integrators[i].param_change(this->T, this->alpha, this->h,
				      this->box_width);
  }
}

void LennardJonesSystem::param_set(float T, float alpha, float h,
				   float box_width, float m, float rmin,
				   float rc, float epsilon) {
  this->box_width = box_width;

  this->m = m;
  this->rmin = rmin;
  this->rc = rc;
  this->epsilon = epsilon;

  this->T = T;
  this->alpha = alpha;
  this->h = h;
}

void LennardJonesSystem::step() {
  // cout << this->nparticles << endl;
  // Set all accelerations to zero
  for (int i = 0; i < this->nparticles; i++) {
    this->accels[i].x = 0.;
    this->accels[i].y = 0.;
  }

  // Calculate accelerations for each point
  for (int i = 0; i < nparticles - 1; i++) {
    for (int j = i + 1; j < nparticles; j++) {
      // TODO: check the sign here
      twofloat ijaccel = acceleration(this->centers[j], this->centers[i]);
      this->accels[i].x += ijaccel.x;
      this->accels[i].y += ijaccel.y;
      // Opposite sign for particle j
      this->accels[j].x -= ijaccel.x;
      this->accels[j].y -= ijaccel.y;
    }
    // Particle i acceleration complete
  }

  // Now update everything
  // The last particle has the correct acceleration after the loop
  for (int i = 0; i < nparticles; i++) {
    this->integrators[i].update(this->accels[i]);
  }
  //this->integrators[nparticles - 1].update(this->accels[nparticles - 1]);
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
  if (-1 * dx > this->box_width / 2) {
    dx += this->box_width;
  }
  if (dy > this->box_width / 2) {
    dy -= this->box_width;
  }
  if (-1 * dy > this->box_width / 2) {
    dy += this->box_width;
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

/**
 * Plain c API for interfacing with python ctypes
 *
 * Some relevant references for doing this:
 * https://stackoverflow.com/questions/1615813/how-to-use-c-classes-with-ctypes/7061012#7061012
 * http://bigbang.waterlin.org/bang/using-python-ctypes-to-link-cpp-library/
 */
extern "C" {
  LennardJonesSystem *LennardJones_new_full(int nparticles, float *centers,
					    float m, float rmin, float epsilon,
					    float T, float alpha, float h,
					    float rc, float box_width) {
    LennardJonesSystem *sys = new LennardJonesSystem(nparticles, centers, m,
						     rmin, epsilon, T,
						     alpha, h, rc, box_width);
    return sys;
  }

  LennardJonesSystem *LennardJones_new(int nparticles, float *centers,
				       float m, float rmin, float epsilon,
				       float T, float alpha, float h,
				       float box_width) {
    LennardJonesSystem *sys = new LennardJonesSystem(nparticles, centers, m,
						     rmin, epsilon, T,
						     alpha, h, box_width);
    return sys;
  }

  float LennardJones_get_box_width(LennardJonesSystem *sys) {
    return sys->get_box_width();
  }

  float LennardJones_get_T(LennardJonesSystem *sys) {
    return sys->get_T();
  }

  float LennardJones_get_alpha(LennardJonesSystem *sys) {
    return sys->get_alpha();
  }

  float LennardJones_get_h(LennardJonesSystem *sys) {
    return sys->get_h();
  }

  float LennardJones_get_m(LennardJonesSystem *sys) {
    return sys->get_m();
  }

  float LennardJones_get_rmin(LennardJonesSystem *sys) {
    return sys->get_rmin();
  }

  float LennardJones_get_rc(LennardJonesSystem *sys) {
    return sys->get_rc();
  }

  float LennardJones_get_epsilon(LennardJonesSystem *sys) {
    return sys->get_epsilon();
  }

  float LennardJones_get_nparticles(LennardJonesSystem *sys) {
    return sys->get_nparticles();
  }

  twofloat *LennardJones_get_centers(LennardJonesSystem *sys) {
    return sys->get_centers();
  }
  
  void LennardJones_param_change_full(LennardJonesSystem *sys, float T,
				      float alpha, float h, float box_width,
				      float m, float rmin, float rc,
				      float epsilon) {
    sys->param_change(T, alpha, h, box_width, m, rmin, rc, epsilon);
  }

  void LennardJones_param_change(LennardJonesSystem *sys, float T, float alpha,
				 float h, float box_width) {
    sys->param_change(T, alpha, h, box_width);
  }

  void LennardJones_step(LennardJonesSystem *sys) {
    sys->step();
  }

  void LennardJones_delete(LennardJonesSystem *sys) {
    delete sys;
  }
}

#endif
