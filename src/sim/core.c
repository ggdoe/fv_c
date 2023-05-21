#include "sim.h"

struct grid grid;

struct pstate pstate;
struct cstate cstate;
struct fluxes fluxes_x;
struct fluxes fluxes_y;

void swap_uv(struct pcell *q);

void update_cells()
{
    size_t l = grid.Nx_tot; // access top/bot cell 
    size_t lf = grid.Nx_tot - 1; // access top/bot cell for fluxes

    // #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 1; j++)
        for (size_t i = 0; i < grid.Nx_tot - 1; i++)
        {
            size_t id = j * grid.Nx_tot + i;
            size_t idf = j * (grid.Nx_tot - 1) + i;

            struct pcell pl;
            struct pcell pr;
            struct fcell F;

            // --- x axis
            pl = get_cell(&pstate, id);
            pr = get_cell(&pstate, id+1);
            // reconstruct_interface_x(&pl, id,   +1.0);
            // reconstruct_interface_x(&pr, id+1, -1.0);

            F = solve_fluxes(&pl, &pr);
            fluxes_x.r[idf]  = F.r;
            fluxes_x.ru[idf] = F.ru;
            fluxes_x.rv[idf] = F.rv;
            fluxes_x.e[idf]  = F.e;
            
            // --- y axis
            pl = get_cell(&pstate, id);
            pr = get_cell(&pstate, id+l);
            // reconstruct_interface_y(&pl, id,   +1.0);
            // reconstruct_interface_y(&pr, id+l, -1.0);
            
            swap_uv(&pl);
            swap_uv(&pr);
            F = solve_fluxes(&pl, &pr);
            fluxes_y.r[idf]  = F.r;
            fluxes_y.ru[idf] = F.rv; // swap
            fluxes_y.rv[idf] = F.ru; // swap
            fluxes_y.e[idf]  = F.e;
        }
        
    // #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 2; j++)
        for (size_t i = 0; i < grid.Nx_tot - 2; i++)
        {
            real dtdx = grid.dt / grid.dx;

            size_t id  = j * grid.Nx_tot + i;
            size_t idf = j * (grid.Nx_tot-1) + i;
           
            cstate.r[id+1]   += dtdx * (fluxes_x.r[idf]  - fluxes_x.r[idf+1]);
            cstate.ru[id+1]  += dtdx * (fluxes_x.ru[idf] - fluxes_x.ru[idf+1]);
            cstate.rv[id+1]  += dtdx * (fluxes_x.rv[idf] - fluxes_x.rv[idf+1]);
            cstate.e[id+1]   += dtdx * (fluxes_x.e[idf]  - fluxes_x.e[idf+1]);
        }

    // #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 2; j++)
        for (size_t i = 0; i < grid.Nx_tot - 2; i++)
        {
            real dtdy = grid.dt / grid.dy;

            size_t id  = j * grid.Nx_tot + i;
            size_t idf = j * (grid.Nx_tot-1) + i;
           
            cstate.r[id+lf]  += dtdy * (fluxes_y.r[idf]  - fluxes_y.r[idf+lf]);
            cstate.ru[id+lf] += dtdy * (fluxes_y.ru[idf] - fluxes_y.ru[idf+lf]);
            cstate.rv[id+lf] += dtdy * (fluxes_y.rv[idf] - fluxes_y.rv[idf+lf]);
            cstate.e[id+lf]  += dtdy * (fluxes_y.e[idf]  - fluxes_y.e[idf+lf]);
        }
}

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

void compute_dt()
{
    real dt = __FLT_MAX__;

    // TODO : attention aux ghost-cells ?
    // #pragma omp parallel for reduction(min:dt)
    for (size_t i = 0; i < grid.N_tot; i++)
    {
        real cs = SQRT(pstate.p[i] * grid.gamma / pstate.r[i]);
        real dtx = grid.dx / (cs + fabs(pstate.u[i]));
        real dty = grid.dy / (cs + fabs(pstate.v[i]));
        real new_dt = (dtx < dty) ? dtx : dty;
        dt = (new_dt < dt) ? new_dt : dt;
    }
    grid.dt = grid.CFL * dt;
    // printf("%.10lf\n", grid.dt);
}

void step(real dt_max)
{
    fill_boundaries();
    conservative_to_primitive(&pstate, &cstate);
    
    compute_dt();
    if(grid.dt > dt_max)
        grid.dt = dt_max;
    grid.t += grid.dt;

    compute_slopes();
    update_cells();
}

void swap_uv(struct pcell *q)
{
    real swap = q->u;
    q->u = q->v;
    q->v = swap;
}
