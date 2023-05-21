#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;
extern struct cstate cstate;
extern struct fluxes fluxes_x;
extern struct fluxes fluxes_y;

extern struct pstate x_slope;
extern struct pstate y_slope;

void linspace(real* dest, real xmin, real xmax, u32 n);
void malloc_struct(void *in, size_t n);

void init_sim(u32 nx, u32 ny)
{
    u32 Ng = 2;
    grid.Ng = Ng;
    grid.CFL = 0.3;
    grid.gamma = 1.4;

    // x discretization
    grid.Nx_tot = nx + 2 * Ng;
    real xmin = 0.; grid.xmin = xmin;
    real xmax = 1.; grid.xmax = xmax;
    real dx = (xmax - xmin) / nx;
    grid.dx = dx;
    grid.x = malloc(grid.Nx_tot * sizeof(real));
    linspace(grid.x, xmin + (0.5 - (real)Ng) * dx, xmax + ((real)Ng - 0.5) * dx, grid.Nx_tot);

    // y discretization
    grid.Ny_tot = ny + 2 * Ng;
    real ymin = 0.; grid.ymin = ymin;
    real ymax = 1.; grid.ymax = ymax;
    real dy = (ymax - ymin) / ny;
    grid.dy = dy;
    grid.y = malloc(grid.Ny_tot * sizeof(real));
    linspace(grid.y, ymin + (0.5 - (real)Ng) * dy, ymax + ((real)Ng - 0.5) * dy, grid.Ny_tot);
    
    grid.N_tot = grid.Nx_tot * grid.Ny_tot;

    malloc_struct(&pstate, grid.Nx_tot * grid.Ny_tot);
    malloc_struct(&cstate, grid.Nx_tot * grid.Ny_tot);
    malloc_struct(&fluxes_x, (grid.Nx_tot - 1) * (grid.Ny_tot - 1));
    malloc_struct(&fluxes_y, (grid.Nx_tot - 1) * (grid.Ny_tot - 1));
    malloc_struct(&x_slope, grid.Nx_tot * grid.Ny_tot);
    malloc_struct(&y_slope, grid.Nx_tot * grid.Ny_tot);

    init_problem(&pstate);

    // memset(fluxes.r, 0, (grid.Nx_tot - 1) * sizeof(real));
    // memset(fluxes.ru, 0, (grid.Nx_tot - 1) * sizeof(real));
    // memset(fluxes.e, 0, (grid.Nx_tot - 1) * sizeof(real));
}

void malloc_struct(void *in, size_t n)
{
    struct pstate *q = in;
    q->r = malloc(n * sizeof(real));
    q->u = malloc(n * sizeof(real));
    q->v = malloc(n * sizeof(real));
    q->p = malloc(n * sizeof(real));
     
    // TODO suppr
    // memset(q->r, 0, n * sizeof(real));
    // memset(q->u, 0, n * sizeof(real));
    // memset(q->v, 0, n * sizeof(real));
    // memset(q->p, 0, n * sizeof(real));
}

void close_sim()
{
    free(grid.x);
    free(grid.y);

    free(pstate.r);
    free(pstate.u);
    free(pstate.v);
    free(pstate.p);

    free(cstate.r);
    free(cstate.ru);
    free(cstate.rv);
    free(cstate.e);

    free(fluxes_x.r);
    free(fluxes_x.ru);
    free(fluxes_x.rv);
    free(fluxes_x.e);

    free(fluxes_y.r);
    free(fluxes_y.ru);
    free(fluxes_y.rv);
    free(fluxes_y.e);

    free(x_slope.r);
    free(x_slope.u);
    free(x_slope.v);
    free(x_slope.p);   

    free(y_slope.r);
    free(y_slope.u);
    free(y_slope.v);
    free(y_slope.p);    
}

