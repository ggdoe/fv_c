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

struct ccell{
    real r;
    real ru;
    real e;
};

struct pcell{
    real r;
    real u;
    real p;
};

struct fcell{
    real r;
    real ru;
    real e;
};


struct pstate{
    real *r;
    real *u;
    real *p;
};

struct cstate{
    real *r;
    real *ru;
    real *e;
};

struct fluxes{
    real *r;
    real *ru;
    real *e;
};

struct grid {
    u32 Ng;
    u32 Nx;
    u32 Nx_tot;
    real dx;
    real *x;

    real xmin, xmax;
    
    real gamma;
    real CFL; // TODO: move
    real dt;
    real t;
};


void init_sim(u32 Nx);
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

// tools
void linspace(real* dest, real xmin, real xmax, u32 n);
void print_mat(real* mat, u32 n, u32 m);


