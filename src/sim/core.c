#include "sim.h"

struct grid grid;

struct pstate pstate;
struct cstate cstate;

void run(real tmax)
{
    real *t = &grid.t;
    // size_t itt = 0;

    while(*t < tmax)
    {
        step(tmax-*t);
        // itt++;
        // printf("%lu - %lf\n", itt, t);
    }

}

void compute_dt(struct pstate *pstate)
{
    real dt = __FLT_MAX__;

    // TODO : attention aux ghost-cells ?
    #pragma omp parallel for reduction(min:dt)
    for (size_t i = 0; i < grid.N_tot; i++)
    {
        real cs = SQRT(pstate->p[i] * grid.gamma / pstate->r[i]);
        real dtx = grid.dx / (cs + fabs(pstate->u[i]));
        real dty = grid.dy / (cs + fabs(pstate->v[i]));
        real new_dt = (dtx < dty) ? dtx : dty;
        dt = (new_dt < dt) ? new_dt : dt;
    }
    grid.dt = grid.CFL * dt;
    // printf("%.10lf\n", grid.dt);
}

void step(real dt_max)
{
    compute_dt(&pstate);
    if(grid.dt > dt_max)
        grid.dt = dt_max;
    grid.t += grid.dt;

    compute_slopes(&pstate);
    update_cells(&pstate, &cstate);
    fill_boundaries(&cstate);
    conservative_to_primitive(&cstate, &pstate);
    
}
