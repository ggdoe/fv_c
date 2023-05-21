#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;
extern struct cstate cstate;
extern struct fluxes fluxes;

extern struct pstate x_slope;

void linspace(real* dest, real xmin, real xmax, u32 n);

void init_sim(u32 Nx)
{
    u32 Ng = 2;

    grid.Ng = Ng;
    grid.Nx = Nx;
    grid.Nx_tot = Nx + 2 * Ng;

    real xmin = 0.; grid.xmin = xmin;
    real xmax = 1.; grid.xmax = xmax;
    real dx = (xmax - xmin) / Nx;
    grid.dx = dx;
    grid.CFL = 0.3;
    grid.gamma = 1.4;

    grid.x = malloc(grid.Nx_tot * sizeof(real));
    linspace(grid.x, xmin + (0.5 - (real)Ng) * dx, xmax + ((real)Ng - 0.5) * dx, grid.Nx_tot);

    pstate.r  = malloc(grid.Nx_tot * sizeof(real));
    pstate.u  = malloc(grid.Nx_tot * sizeof(real));
    pstate.p  = malloc(grid.Nx_tot * sizeof(real));
    cstate.r  = malloc(grid.Nx_tot * sizeof(real));
    cstate.ru = malloc(grid.Nx_tot * sizeof(real));
    cstate.e  = malloc(grid.Nx_tot * sizeof(real));
    fluxes.r  = malloc((grid.Nx_tot - 1) * sizeof(real));
    fluxes.ru = malloc((grid.Nx_tot - 1) * sizeof(real));
    fluxes.e  = malloc((grid.Nx_tot - 1) * sizeof(real));

    x_slope.r = malloc(grid.Nx_tot * sizeof(real));
    x_slope.u = malloc(grid.Nx_tot * sizeof(real));
    x_slope.p = malloc(grid.Nx_tot * sizeof(real));

    init_problem(&pstate);

    // memset(fluxes.r, 0, (grid.Nx_tot - 1) * sizeof(real));
    // memset(fluxes.ru, 0, (grid.Nx_tot - 1) * sizeof(real));
    // memset(fluxes.e, 0, (grid.Nx_tot - 1) * sizeof(real));
}

void close_sim()
{
    free(grid.x);
    free(pstate.r);
    free(pstate.u);
    free(pstate.p);

    free(cstate.r);
    free(cstate.ru);
    free(cstate.e);

    free(fluxes.r);
    free(fluxes.ru);
    free(fluxes.e);

    free(x_slope.r);
    free(x_slope.u);
    free(x_slope.p);    
}

