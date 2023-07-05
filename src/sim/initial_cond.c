#include "sim.h"

#ifdef _WIN32
#define __attribute_maybe_unused__
#endif

extern struct grid grid;

__attribute_maybe_unused__ 
void init_SOD_2d(struct pstate *state)
{
    size_t Nx = grid.Nx_tot;
    size_t Ny = grid.Ny_tot;
    
    // SOD-2d
    for (size_t j = 0; j < Ny; j++)
    {
        for (size_t i = 0; i < Nx / 2; i++)
        {
            state->r[i * Nx + j] = 1.0;
            state->u[i * Nx + j] = 0.;
            state->v[i * Nx + j] = 0.;
            state->p[i * Nx + j] = 1.0;
        }

        for (size_t i = Nx / 2; i < Nx; i++)
        {
            state->r[i * Nx + j] = 0.125;
            state->u[i * Nx + j] = 0.;
            state->v[i * Nx + j] = 0.;
            state->p[i * Nx + j] = 0.1;
        }
    }
}

__attribute_maybe_unused__ 
void init_kelvin_helmholtz(struct pstate *pstate)
{
    size_t Nx = grid.Nx_tot;
    size_t Ny = grid.Ny_tot;

    srand(time(NULL));
    
    // SOD-2d
    for (size_t j = 0; j < Ny; j++)
        for (size_t i = 0; i < Nx; i++)
        {
            real r,u;
            if(j <= (Ny+2) / 4 || j > 3 * (Ny-2) / 4){
                r = 2.0; u = -0.5;
            }
            else{
                r = 1.0; u = 0.5;
            }

            pstate->r[j * Nx + i] = r;
            pstate->u[j * Nx + i] = u;
            pstate->v[j * Nx + i] = 0.0;
            pstate->p[j * Nx + i] = 2.5;

            pstate->u[j * Nx + i] += 0.01 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
            pstate->v[j * Nx + i] += 0.01 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
            // if(i < Nx / 2)
            //     state->u[j * Nx + i] += +0.01;
            // else
            //     state->v[j * Nx + i] += -0.01;
        }
}

__attribute_maybe_unused__ 
void init_blast(struct pstate *pstate)
{
    size_t Nx = grid.Nx_tot;
    size_t Ny = grid.Ny_tot;

    real x_shift = 0.5;
    real y_shift = 0.5;
    real radius  = 0.2;
    real p_in    = 10.0;
    real p_out   = 1.0;
    real rho_in  = 1.0;
    real rho_out = 10.0;

    const real thck = 0.02; // blast smoothing
    
    // SOD-2d
    for (size_t j = 0; j < Ny; j++)
        for (size_t i = 0; i < Nx; i++)
        {
            const real x = grid.x[i] - x_shift;
            const real y = grid.y[j] - y_shift;
            const real r = sqrt(x*x + y*y);

            real tr = 0.5 * (tanh((r - radius) / thck) + 1.0);
            // tr = (r < radius) ? 0.0 : 1.0;

            pstate->r[j * Nx + i] = rho_out * tr + rho_in * (1-tr);
            pstate->u[j * Nx + i] = 0; // 0.001 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
            pstate->v[j * Nx + i] = 0; // 0.001 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
            pstate->p[j * Nx + i] = p_out * tr + p_in * (1-tr);
        }
}

void init_problem(struct pstate *pstate, struct cstate *cstate)
{
    grid.t = 0;

    init_blast(pstate);
    // init_SOD_2d(pstate);
    // init_kelvin_helmholtz(pstate);
    primitive_to_conservative(pstate, cstate);
    fill_boundaries(cstate);
}

