/****************************************************************************
 * yang-mills.cpp
 * David Boozer
 * 14 September 2009
 * 26 March 2011 (convert to c++)
 * 26 April 2013 (see old-yang-mills.cpp)
 *****************************************************************************
 * Simulate SU(2) gauge theory in (1+1) dimensions in Lorentz gauge.
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <string>
#include <stdio.h>
#include <time.h>

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

const int N = 2001;                          // number of gridpoints
const int NUM_PARTICLES = 4;                  // number of particles
const double L = 2.0;                        // size of region
const double G = 23.75;                        // gauge field coupling
const double DELTA = L/(N-1);                // spacing of gridpoints
const double delta = 0.01;                  // spacing of time sampling points
// const double SIGMA = 0.001;                  // size of particle
const double SIGMA = DELTA;
const double STEP_MAX = 0.001;               // maximum allowed timestep
const double tmax = 100.0;
const double mass = 0.510999;
const double z0 = 0.02;
const double u0 = 0.1;
const double q0 = 0.303;

// charge density profile of particle
// SIGMA is size of particle, here equal to dx
#define GAUSSIAN(x) (exp(-(x)*(x)/(2*SIGMA*SIGMA))/(SIGMA*sqrt(2*M_PI)))

// DIRAC DELTA FUNCTION
// fabs() is the absolute value function
// if |x|<particle_size then D(x)=gaussian(x), otherwise D(x)=0
#define D(x) ((fabs((x)/SIGMA) < 5.0) ? GAUSSIAN((x)) : 0.0)

// this is an array
double x[N];   // location of gridpoints

class vector {
public:
  double x;
  double y;
  double z;

  void set_to_zero ();
  vector operator- ();
  vector& operator+= (vector);
};

// state of whole system i.e. all particles across all space grid
class state {
public:
  // remember, N is the amount of grid points so these are N-dimensional lists corresponding to x-length
  vector A[N], V[N], E[N], q[NUM_PARTICLES];
  double z[NUM_PARTICLES], u[NUM_PARTICLES];

  state dot ();

  void evolve (double);
  void initialize ();
  void gauss_violation (vector*);

  double energy_particle ();
  double energy_field ();
  double energy_total ();
  double beta (int);

  vector rho (int);
  vector J (int);
  vector field_on_wl (vector*, int);

  void output_frame ();
};

/****************************************************************************
 * vector operations
 ****************************************************************************/
// HOW DOES THE VECTOR CLASS KNOW TO INCLUDE THESE FUNCTIONS?

// vector inversion (parity operator on the vector)
vector vector::operator- () {
  vector b;

  // this is a function inside of the vector class ("::") so it already knows what x,y,z are
  b.x = -x;
  b.y = -y;
  b.z = -z;

  return b;
}

// vector cross product
vector operator^ (vector a, vector b) {
  vector c;

  c.x = a.y*b.z - a.z*b.y;
  c.y = a.z*b.x - a.x*b.z;
  c.z = a.x*b.y - a.y*b.x;

  return c;
}

// vector addition
vector operator+ (vector a, vector b) {
  vector c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;

  return c;
}

// vector subtraction
vector operator- (vector a, vector b) {
  vector c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;

  return c;
}

// vector dot product
double operator* (vector a, vector b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

// vector multiplication by scalar
vector operator* (double a, vector b) {
  vector c;

  c.x = a*b.x;
  c.y = a*b.y;
  c.z = a*b.z;

  return c;
}

// vector multiplication by scalar (other way round)
vector operator* (vector a, double b) {
  vector c;

  c.x = a.x*b;
  c.y = a.y*b;
  c.z = a.z*b;

  return c;
}

// vector division by scalar from right
vector operator/ (vector a, double b) {
  vector c;

  c.x = a.x/b;
  c.y = a.y/b;
  c.z = a.z/b;

  return c;
}

// vector += operation
vector& vector::operator+= (vector a) {
  x += a.x;
  y += a.y;
  z += a.z;

  // what is the *this pointer?
  return *this;
}

void vector::set_to_zero () {
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

/****************************************************************************
 * operator overloading
 ****************************************************************************/

// addition of two states
state operator+ (state s1, state s2) {
  state s_sum;
  int n, i;

  for (n=0; n<N; n++) {
    s_sum.A[n] = s1.A[n] + s2.A[n];
    s_sum.V[n] = s1.V[n] + s2.V[n];
    s_sum.E[n] = s1.E[n] + s2.E[n];
  }

  for (i=0; i<NUM_PARTICLES; i++) {
    s_sum.z[i] = s1.z[i] + s2.z[i];
    s_sum.u[i] = s1.u[i] + s2.u[i];
    s_sum.q[i] = s1.q[i] + s2.q[i];
  }

  return s_sum;
}

// multiplies all variables by a scalar c
state operator* (double c, state s) {
  state s_product;
  int n, i;

  for (n=0; n<N; n++) {
    s_product.A[n] = c*s.A[n];
    s_product.V[n] = c*s.V[n];
    s_product.E[n] = c*s.E[n];
  }

  for (i=0; i<NUM_PARTICLES; i++) {
    s_product.z[i] = c*s.z[i];
    s_product.u[i] = c*s.u[i];
    s_product.q[i] = c*s.q[i];
  }

  return s_product;
}

/****************************************************************************
 * lattice derivative
 ****************************************************************************/

// differentiates a and returns into a_x
void diff_x (vector a[N], vector a_x[N]) {
  int n;

  // central difference
  for (n=1; n<N-1; n++) a_x[n] = (a[n+1] - a[n-1])/(2*DELTA);

  // boundary conditions
  a_x[0] = (a[1] - a[0])/DELTA;     // foward difference???
  a_x[N-1] = (a[N-1] - a[N-2])/DELTA;
}

/****************************************************************************
 * manipulate states
 ****************************************************************************/

// propagate in time
// tau is the timestep to evolve for
void state::evolve (double tau) {
  state dot1, dot2, dot3, dot4;
  int n, num_steps;
  double step;

  num_steps = (int) ceil (fabs(tau)/STEP_MAX);
  step = tau/num_steps;

  // loop across time points inside tau
  // the "*this" pointer is the same as the "self" in Python classes
  // Why does he do time propagation like this?
    // i think it is just a way to approximate the time interval by "averaging"
    // the change over the interval (so the change is not dS at left point of grid step
    // but actually an average of dS over all of the time step)
  for (n=0; n<num_steps; n++) {
    dot1 = (*this).dot();
    dot2 = (*this + 0.5*step*dot1).dot();
    dot3 = (*this + 0.5*step*dot2).dot();
    dot4 = (*this + step*dot3).dot();
    *this = *this + (step/6.0)*(dot1 + 2*dot2 + 2*dot3 + dot4);
  }
}

// just sets all variables of an instance of the state class to zero
void state::initialize () {
  int n, i;

  for (i=0; i<NUM_PARTICLES; i++) {
    z[i] = 0.0;
    u[i] = 0.0;
    q[i].set_to_zero();
  }

  for (n=0; n<N; n++) {
    A[n].set_to_zero();
    V[n].set_to_zero();
    E[n].set_to_zero();
  }
}

/****************************************************************************
 * properties of states
 ****************************************************************************/

// evaluates classical velocity of particle i
double state::beta (int i) {
  return u[i]/sqrt(1 + u[i]*u[i]);
}

// evaluates the charge density at grid position x[n]
vector state::rho (int n) {
  vector temp;
  int i;

  temp.set_to_zero();
  for (i=0; i<NUM_PARTICLES; i++) temp += q[i]*D(x[n] - z[i]);
  return temp;
}


vector state::J (int n) {
  vector temp;
  int i;

  temp.set_to_zero();
  for (i=0; i<NUM_PARTICLES; i++) temp += q[i]*beta(i)*D(x[n] - z[i]);
  return temp;
}

// field on world line?
vector state::field_on_wl (vector field[N], int i) {
  vector temp;
  int n;

  temp.set_to_zero();
  // loop across grid points
  for (n=0; n<N; n++) temp += field[n]*D(x[n] - z[i])*DELTA;
  return temp;
}

double state::energy_particle () {
  double energy;
  int i;

  energy = 0.0;
  for (i=0; i<NUM_PARTICLES; i++) energy += sqrt(1 + u[i]*u[i]);
  return energy - NUM_PARTICLES;
}

double state::energy_field () {
  double energy;
  int n;

  energy = 0.0;
  for (n=0; n<N; n++) energy += 0.25*E[n]*E[n]*DELTA;
  return energy;
}

double state::energy_total () {
  return energy_field() + energy_particle();
}

// calculate the constraint vector g
void state::gauss_violation (vector g[N]) {
  vector E_x[N];
  int n;

  diff_x (E, E_x);
  // loop over spatial grid
  for (n=0; n<N; n++) g[n] = 2*rho(n) + G*(A[n]^E[n]) - E_x[n];
}

/****************************************************************************
 * equations of motion
 ****************************************************************************/

#define ALPHA 0.5
#define BETA 0.5

state state::dot () {
  state dot;
  vector V_x[N], A_x[N];
  int n, i;

  diff_x (V, V_x);
  diff_x (A, A_x);

  for (i=0; i<NUM_PARTICLES; i++) {
    dot.u[i] = q[i]*field_on_wl(E,i)/mass;
    dot.z[i] = beta(i);
    dot.q[i] = -G*(field_on_wl(V,i) - beta(i)*field_on_wl(A,i))^q[i];
  }

  for (n=1; n<N-1; n++) {
    dot.E[n] = -2*J(n) - G*(V[n]^E[n]);
    dot.A[n] = -V_x[n] - E[n] + G*(A[n]^V[n]);
    dot.V[n] = -A_x[n];
  }

  // BOUNDARY CONDITIONS !
  dot.E[0].set_to_zero();
  dot.E[N-1].set_to_zero();

  dot.A[0] = -ALPHA*V_x[0] + BETA*A_x[0];
  dot.A[N-1] = -ALPHA*V_x[N-1] - BETA*A_x[N-1];

  dot.V[0] = -ALPHA*A_x[0] + BETA*V_x[0];
  dot.V[N-1] = -ALPHA*A_x[N-1] - BETA*V_x[N-1];

  return dot;
}

/****************************************************************************
 * experiments
 ****************************************************************************/

void state::output_frame () {
  FILE *file_field;
  char filename[32];
  static int frame_num = 0;
  int n;

  sprintf (filename, "frame-%04d.txt", frame_num++);
  file_field = fopen (filename, "w");
  for (n=0; n<N; n++)
    fprintf (file_field, "%f, %f, %f, %f\n",
	     x[n], E[n].x, V[n].x, A[n].x);
  fclose (file_field);
}

void test () {
  FILE *file_particle;
  FILE *file_params;
  state state;
  double t;
  int n;

  state.initialize();

  state.z[0] = z0;
  state.z[1] = z0;
  state.z[2] = -z0;
  state.z[3] = -z0;

  state.q[0].x = -q0;
  state.q[1].x = q0;
  state.q[2].y = -q0;
  state.q[3].y = q0;

  state.u[0] = u0;
  state.u[1] = -u0;
  state.u[2] = u0;
  state.u[3] = -u0;

/*
  state.output_frame ();
  file_particle = fopen ("particle.txt", "w");
  // loop over times
  for (t=0; t<10.0; t += delta) {
    // save data
    fprintf (file_particle, "%f, %f, %f, %f, %f, %f, %f\n", t,
	     state.u[0], state.z[0], state.q[0].x,
	     state.energy_particle (), state.energy_field (),
	     state.energy_total ());
*/
  //char* filename = concat(concat("data",time(0)),".txt");
  file_particle = fopen ("data.txt", "w");
  file_params = fopen ("dataparams.txt", "w");

  // save parametres to "dataparams.txt"
  fprintf (file_params, "%f\n", G);


  // loop over times
  for (t=0; t<tmax; t+= delta) {
    // save data as t, z1, z2, z3, z4, u1, u2, u3, u4, q11, q12, q13, q21, q22, q23, q31, q32, q33, q41, q42, q43,
    fprintf (file_particle, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",
                            t, state.z[0], state.z[1], state.z[2], state.z[3],
                               state.u[0], state.u[1], state.u[2], state.u[3],
                               state.q[0].x, state.q[0].y, state.q[0].z,
                               state.q[1].x, state.q[1].y, state.q[1].z,
                               state.q[2].x, state.q[2].y, state.q[2].z,
                               state.q[3].x, state.q[3].y, state.q[3].z);


    // evolve dt
    state.evolve (delta);
  }
  fclose (file_particle);
  //state.output_frame ();
}

/****************************************************************************
 * generate output files
 ****************************************************************************/

int main () {
  int n;

  // create spatial vector
  for (n=0; n<N; n++) x[n] = 0.5*L*((double) 2*n/(N-1) - 1);

  // execute the function
  test ();
}
