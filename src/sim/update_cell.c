#include "sim.h"

extern struct grid grid;
struct fluxes fluxes_x;
struct fluxes fluxes_y;

static void swap_uv(struct pcell *q)
{
    real swap = q->u;
    q->u = q->v;
    q->v = swap;
}

// TODO : 
//        - passer en full SOA

void update_cells(struct pstate *pstate, struct cstate *cstate)
{
    size_t l = grid.Nx_tot; // access top/bot cell 

    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot; j++)
        for (size_t i = 0; i < grid.Nx_tot - 1; i++)
        {
            size_t id = j * l + i;
            struct pcell pl;
            struct pcell pr;
            struct fcell F;

            // --- x axis
            pl = get_cell(pstate, id);
            pr = get_cell(pstate, id+1);
            reconstruct_muscl_hancock(&pl, id);
            reconstruct_muscl_hancock(&pr, id+1);
            reconstruct_interface_x(&pl, id,   +1.0);
            reconstruct_interface_x(&pr, id+1, -1.0);


            F = solve_fluxes(&pl, &pr);
            fluxes_x.r [id] = F.r;
            fluxes_x.ru[id] = F.ru;
            fluxes_x.rv[id] = F.rv;
            fluxes_x.e [id] = F.e;
        }
            
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 1; j++)
        for (size_t i = 0; i < grid.Nx_tot; i++)
        {
            size_t id = j * l + i;
            struct pcell pl;
            struct pcell pr;
            struct fcell F;

            // --- y axis
            pl = get_cell(pstate, id);
            pr = get_cell(pstate, id+l);
            reconstruct_muscl_hancock(&pl, id);
            reconstruct_muscl_hancock(&pr, id+l);
            reconstruct_interface_y(&pl, id,   +1.0);
            reconstruct_interface_y(&pr, id+l, -1.0);
            
            swap_uv(&pl);
            swap_uv(&pr);
            F = solve_fluxes(&pl, &pr);
            fluxes_y.r [id]  = F.r;
            fluxes_y.ru[id] = F.rv; // swap
            fluxes_y.rv[id] = F.ru; // swap
            fluxes_y.e [id]  = F.e;
        }
        
    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 2; j++)
        for (size_t i = 0; i < grid.Nx_tot - 2; i++)
        {
            real dtdx = grid.dt / grid.dx;

            size_t id  = j * grid.Nx_tot + i;
           
            cstate->r [id+1]  += dtdx * (fluxes_x.r [id] - fluxes_x.r [id+1]);
            cstate->ru[id+1]  += dtdx * (fluxes_x.ru[id] - fluxes_x.ru[id+1]);
            cstate->rv[id+1]  += dtdx * (fluxes_x.rv[id] - fluxes_x.rv[id+1]);
            cstate->e [id+1]  += dtdx * (fluxes_x.e [id] - fluxes_x.e [id+1]);
        }

    #pragma omp parallel for collapse(2)
    for (size_t j = 0; j < grid.Ny_tot - 2; j++)
        for (size_t i = 0; i < grid.Nx_tot - 2; i++)
        {
            real dtdy = grid.dt / grid.dy;

            size_t id  = j * grid.Nx_tot + i;
           
            cstate->r [id+l] += dtdy * (fluxes_y.r [id] - fluxes_y.r [id+l]);
            cstate->ru[id+l] += dtdy * (fluxes_y.ru[id] - fluxes_y.ru[id+l]);
            cstate->rv[id+l] += dtdy * (fluxes_y.rv[id] - fluxes_y.rv[id+l]);
            cstate->e [id+l] += dtdy * (fluxes_y.e [id] - fluxes_y.e [id+l]);
        }
}
