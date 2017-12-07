/**
 * 
 */
#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

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
  Langevin(float T, float alpha, float h, Point *p);
  ~Langevin();
  void param_change(float T, float alpha, float h);
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
  assert(alpha > 0);
  assert(h > 0);

  param_change(T, alpha, h);

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
void Langevin::param_change(float T, float alpha, float h) {
  this->T = T;
  this->alpha = alpha;
  this->h = h;
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
		     float rc, float box_width);
  LennardJonesSystem(int nparticles, Point *centers, float m, float rmin,
		     float epsilon, float T, float alpha, float h,
		     float box_width);
  ~LennardJonesSystem();
  twofloat acceleration(twofloat from, twofloat on);
  void param_change(float T, float alpha, float h, float box_width,
		    float m, float rmin, float rc, float epsilon);
  void param_change(float T, float alpha, float h);
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
private:
  float box_width;
  float T, alpha, h;
  float m, rmin, rc, epsilon;

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
  construct_lennard_jones(nparticles, centers, m, rmin, epsilon,
			  T, alpha, h, DEFAULT_RC_MULT * rmin, box_width);
}

LennardJonesSystem::LennardJonesSystem(int nparticles, Point *centers, float m,
				       float rmin, float epsilon,
				       float T, float alpha, float h,
				       float rc, float box_width) {
  construct_lennard_jones(nparticles, centers, m, rmin, epsilon,
			  T, alpha, h, rc, box_width);
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
				    &centers[i]);
  }
  //cout << this->nparticles << endl;
}

LennardJonesSystem::~LennardJonesSystem() {
  delete[] this->accels;
  delete[] this->integrators;
}

void LennardJonesSystem::param_change(float T, float alpha, float h) {
  this->T = T;
  this->alpha = alpha;
  this->h = h;

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
    this->integrators[i].param_change(this->T, this->alpha, this->h);
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
  cout << this->nparticles << endl;
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
    this->integrators[i].update(this->accels[i]);
  }

  // The last particle has the correct acceleration after the loop
  this->integrators[nparticles - 1].update(this->accels[nparticles - 1]);
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

/**
 * Plain c API for interfacing with python ctypes
 *
 * Some relevant references for doing this:
 * https://stackoverflow.com/questions/1615813/how-to-use-c-classes-with-ctypes/7061012#7061012
 * http://bigbang.waterlin.org/bang/using-python-ctypes-to-link-cpp-library/
 */
extern "C" {
  LennardJonesSystem *LennardJones_new_full(int nparticles, Point *centers,
					    float m, float rmin, float epsilon,
					    float T, float alpha, float h,
					    float rc, float box_width) {
    LennardJonesSystem *sys = new LennardJonesSystem(nparticles, centers, m,
						     rmin, epsilon, T,
						     alpha, h, rc, box_width);
    return sys;
  }

  LennardJonesSystem *LennardJones_new(int nparticles, Point *centers,
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
  
  void LennardJones_param_change_full(LennardJonesSystem *sys, float T,
				      float alpha, float h, float box_width,
				      float m, float rmin, float rc,
				      float epsilon) {
    sys->param_change(T, alpha, h, box_width, m, rmin, rc, epsilon);
  }

  void LennardJones_param_change(LennardJonesSystem *sys, float T, float alpha,
				 float h) {
    sys->param_change(T, alpha, h);
  }

  void LennardJones_step(LennardJonesSystem *sys) {
    sys->step();
  }

  void LennardJones_delete(LennardJonesSystem *sys) {
    delete sys;
  }
}

#endif
