#include "sim.h"

#ifdef _WIN32
#define __attribute_maybe_unused__
#endif

extern struct grid grid;
extern struct pstate pstate;
extern struct cstate cstate;

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
            state->r[j * Nx + i] = 1.0;
            state->u[j * Nx + i] = 0.;
            state->v[j * Nx + i] = 0.;
            state->p[j * Nx + i] = 1.0;
        }

        for (size_t i = Nx / 2; i < Nx; i++)
        {
            state->r[j * Nx + i] = 0.125;
            state->u[j * Nx + i] = 0.;
            state->v[j * Nx + i] = 0.;
            state->p[j * Nx + i] = 0.1;
        }
    }
}

__attribute_maybe_unused__ 
void init_kelvin_helmholtz(struct pstate *state)
{
    size_t Nx = grid.Nx_tot;
    size_t Ny = grid.Ny_tot;
    
    // SOD-2d
    for (size_t j = 0; j < Ny; j++)
        for (size_t i = 0; i < Nx; i++)
        {
            real r,u;
            if(j > Ny / 4 && j < 3 * Ny / 4){
                r = 1.0; u = 0.5;
            }
            else{
                r = 2.0; u = -0.5;
            }

            state->r[j * Nx + i] = r;
            state->u[j * Nx + i] = u;
            state->v[j * Nx + i] = 0.0;
            state->p[j * Nx + i] = 2.5;

            state->u[j * Nx + i] += 0.01 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
            state->v[j * Nx + i] += 0.01 * (2.*(real)rand()/(real)RAND_MAX - 1.0);
        }
}

void init_problem(struct pstate *state)
{
    grid.t = 0;

    // init_SOD_2d(state);
    init_kelvin_helmholtz(state);

    fill_boundaries();
    primitive_to_conservative(&cstate, &pstate);
}

