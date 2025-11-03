#ifndef SODE_H
#define SODE_H
#define SODE_V 230529 // version (yy.mm.dd)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// solve non-autonomous problem
// dxdt = f_aux(x,t), where aux are system parameters

typedef enum { SODE_RK56_FB,
               SODE_RK56_CK,
               SODE_RK78_DP } sode_type;
typedef enum { SODE_OK,
               SODE_SUCCESS_TIME,
               SODE_SUCCESS_MONITOR,
               SODE_CONTINUE_GOOD_STEP,
               SODE_CONTINUE_BAD_STEP,
               SODE_FAILED,
               SODE_BAD_FUNC } sode_status;

class sode {
  public:
    // class constructor
    sode(sode_type type, int dim);
    // class destructor
    ~sode();
    // configure stepsizes
    void configure(double h_init, double h_min, double h_max);
    // configure tolerances
    void configure(double tol_sol, double tol_mon, double tol_end, double damp);
    // define dynamical system
    void set_system(int (*system)(double *f, double *x, double t, void *aux));
    // function to detect user defined changes in the orbit
    void set_monitor(bool (*monitor)(double *x, double t, void *aux));
    // set auxiliary pointer for the monitor
    void set_aux(void *aux);
    // run before computing new orbits
    void reset();
    // evolve iteratively system up to t_end or up to a change in monitor
    // (-1: true->false, 0: don't monitor, 1:false->true)
    int evolve(double *x, double *t, double t_end, int chg_mon, void *aux);
    // set verbosity
    void set_verb();

  private:
    int  dim; // the dimension of the dynamical system
    bool verb;

    double h_init;  // initial stepsize when begining orbit
    double h_min;   // the minimum allowed step-size
    double h_max;   // the maximum allowed step-size
    double tol_end; // end time tolerance
    double tol_sol; // the maximum distance between solutions
    double tol_mon; // tolerance on the location of monitor change
    double damp;    // stepsize damping

    int (*system)(double *f, double *x, double t, void *aux) = nullptr;
    bool (*monitor)(double *x, double t, void *aux) = nullptr;
    void *monitor_aux = nullptr;

    sode_type type;      // integration method
    int       q;         // consistency order, err is O(h^p)
    int       s;         // the number of steps in the procedure
    double    h;         // current step-size
    bool      bisecting; // true if currently bisecting event
    double  **k;         // the k's of the integration routine
    double  **a;         // the a's in the Butcher table
    double   *c;         // the c's in the Butcher table
    double   *b1;        // the higher order b's
    double   *b2;        // the lower order b's

    double *x0; // auxiliary vector
    double *x1; // higher order solution
    double *x2; // low order solution

    // calculate the k's for the runge-kutta method
    int calc_ks(double *x, double t, void *aux);
    // combine two vectors
    void vector_combine(double *x_comb, double a, double b,
                        double *xa, double *xb, int size);
    // copy vector x to x_copy
    void vector_copy(double *x_copy, double *x, int size);
    // calculate the Chebyshev distance between two vectors
    double vector_dist_chbyv(double *xa, double *xb, int size);
    // allocate Butcher table for Fehlberg 5(6) method
    void alloc_56FB();
    // allocate Butcher table for Cash-Karp 5(6) method
    void alloc_56CK();
    // allocate Butcher table for Dormand-Prince 7(8) method
    void alloc_RK78();
};

#endif // SODE_H
