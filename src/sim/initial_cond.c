#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;
extern struct cstate cstate;

void init_problem(struct pstate *state)
{
    grid.t = 0;
    
    // SOD
    for (size_t i = 0; i < grid.Nx_tot / 2; i++)
    {
        state->r[i] = 1.0;
        state->u[i] = 0.;
        state->p[i] = 1.0;
    }

    for (size_t i = grid.Nx_tot / 2; i < grid.Nx_tot; i++)
    {
        state->r[i] = 0.125;
        state->u[i] = 0.;
        state->p[i] = 0.1;
    }

    fill_boundaries();
    primitive_to_conservative(&cstate, &pstate);
}

