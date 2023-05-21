#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>
#include <stdbool.h>

#define SQRT(x) sqrt(x)
#define POW(x,y) pow(x,y)
typedef double real;
typedef unsigned int u32;

// conservative cell
struct ccell{
    real r;
    real ru;
    real rv;
    real e;
};

// primitive cell
struct pcell{
    real r;
    real u;
    real v;
    real p;
};

// flux cell
struct fcell{
    real r;
    real ru;
    real rv;
    real e;
};

// primitive grid
struct pstate{
    real *r;
    real *u;
    real *v;
    real *p;
};

// conservative grid
struct cstate{
    real *r;
    real *ru;
    real *rv;
    real *e;
};

// fluxes grid
struct fluxes{
    real *r;
    real *ru;
    real *rv;
    real *e;
};

// grid params
struct grid {
    u32 Ng;

    u32 N_tot;

    u32 Nx_tot;
    real dx;
    real *x;

    u32 Ny_tot;
    real dy;
    real *y;

    real xmin, xmax;
    real ymin, ymax;
    
    real gamma;
    real CFL;
    real dt;
    real t;
};


void init_sim(u32 nx, u32 ny);
void close_sim();

void init_problem(struct pstate *state);

struct fcell solve_fluxes(struct pcell *pl, struct pcell *pr);
void primitive_to_conservative(struct cstate *c, struct pstate *p);
void conservative_to_primitive(struct pstate *p, struct cstate *c);
struct pcell get_cell(struct pstate *pstate, size_t index);
struct fcell get_flux(struct pcell *p);

void step(real dt_max);
void run(real tmax);

void fill_boundaries();
void compute_slopes();
void reconstruct_interface_x(struct pcell *p, size_t i, real sign);
void reconstruct_interface_y(struct pcell *p, size_t i, real sign);

// tools
void print_mat(real* mat, u32 n, u32 m);


