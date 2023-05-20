#include "sim.h"

extern struct grid grid;
extern struct pstate pstate;

struct pcell get_cell(struct pstate *pstate, size_t index)
{
    struct pcell cell = {.r = pstate->r[index], 
                         .u = pstate->u[index],
                         .p = pstate->p[index]};
    return cell;
}


void primitive_to_conservative(struct cstate *c, struct pstate *p)
{
    for (size_t i = 0; i < grid.Nx_tot; i++)
    {
        c->r[i]  = p->r[i];
        c->ru[i] = p->r[i] * p->u[i];
        real Ek  = 0.5 * p->r[i] * (p->u[i] * p->u[i]);
        c->e[i]  = Ek + p->p[i] / (grid.gamma - 1.0);
    }
}

void conservative_to_primitive(struct pstate *p, struct cstate *c)
{
    for (size_t i = 0; i < grid.Nx_tot; i++)
    {
        p->r[i] = c->r[i];
        p->u[i] = c->ru[i] / c->r[i];
        real Ek = 0.5 * (c->ru[i] * c->ru[i]) / c->r[i];
        p->p[i] = (c->e[i] - Ek) * (grid.gamma - 1.0);

        if(p->r[i] < 0.0)
            p->r[i] = 1.0e-5;
        if(p->p[i] < 0.0)
            p->p[i] = 1.0e-5;
        // TODO : call 'primitive_to_conservative' maybe
    }
}

struct fcell get_flux(struct pcell *p)
{
    struct fcell c;
    c.r  = p->r * p->u;
    c.ru = p->r * p->u * p->u + p->p;
    real E = 0.5 * p->r * (p->u * p->u) + p->p / (grid.gamma - 1.0);
    c.e = (E + p->p) * p->u;
    
    return c;
}

