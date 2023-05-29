#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <assert.h>
#include <stdbool.h>

#include "time.h"

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
    u32 Ny_tot;

    real xmin, xmax;
    real ymin, ymax;

    real *x;
    real *y;

    real dx;
    real dy;
    
    real gamma;
    real CFL;
    real dt;
    real t;
};


void init_sim(u32 nx, u32 ny);
void reset_sim();
void close_sim();

void init_problem(struct pstate *pstate, struct cstate *cstate);

struct fcell solve_fluxes(struct pcell *pl, struct pcell *pr);
void primitive_to_conservative(struct pstate *p, struct cstate *c);
void conservative_to_primitive(struct cstate *c, struct pstate *p);
struct pcell get_cell(struct pstate *pstate, size_t index);
struct fcell get_flux(struct pcell *p);

void update_cells(struct pstate *pstate, struct cstate *cstate);
void step(real dt_max);
void run(real tmax);

void fill_boundaries(struct cstate *cstate);
void compute_slopes(struct pstate *pstate);
void reconstruct_interface_x(struct pcell *p, size_t i, real sign);
void reconstruct_interface_y(struct pcell *p, size_t i, real sign);
void reconstruct_muscl_hancock(struct pcell *p, size_t i);

// tools
void print_mat(real* mat, u32 n, u32 m);


