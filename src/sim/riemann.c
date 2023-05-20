#include "sim.h"

extern struct grid grid;

struct interface_val
{
    real L;
    real R;
};

struct ccell cell_primitive_to_conservative(struct pcell *p)
{
    struct ccell c;
    c.r  = p->r;
    c.ru = p->r * p->u;
    real Ek = 0.5 * p->r * (p->u * p->u);
    c.e = Ek + p->p / (grid.gamma - 1.0);

    return c;
}

struct interface_val compute_wave_speed(struct pcell *pl, struct pcell *pr)
{
    // toro2019
    real gamma = grid.gamma;

    real z  = (gamma - 1.0) / (2.0 * gamma);
    
    real aL = SQRT(pl->p * gamma / pl->r);
    real aR = SQRT(pr->p * gamma / pr->r);

    real numerator = (aL + aR - gamma * z * (pr->u - pl->u));
    real denominator = aL / POW(pl->p, z) + aR / POW(pr->p, z);
    // if(numerator <= 1e-2) 
    //     numerator = 1e-2;
    real ps = POW(numerator / denominator, 1.0/z);


    z = (gamma + 1.0) / (2.0 * gamma);
    real qL = (ps < pl->p) ? 1.0 : SQRT(1.0 + z * (ps / pl->p - 1.0));
    real qR = (ps < pr->p) ? 1.0 : SQRT(1.0 + z * (ps / pr->p - 1.0));

    struct interface_val S = {.L = pl->u - aL * qL,
                              .R = pr->u + aR * qR};
    return S;
}

struct fcell get_flux_star_shift(struct ccell *U, struct pcell *Q, real S, real Sstar)
{
    real factor = Q->r * (S - Q->u) / (S - Sstar);

    // Ustar
    struct fcell shift = {.r  = factor * 1.0,
                          .ru = factor * Sstar,
                          .e  = factor * (U->e / U->r 
                                          + (Sstar - Q->u) * (Sstar + Q->p / (Q->r * (S - Q->u))))
                         };
    
    // shift flux
    shift.r  = S * (shift.r  - U->r);
    shift.ru = S * (shift.ru - U->ru);
    shift.e  = S * (shift.e  - U->e);

    return shift;
}

struct fcell solve_fluxes(struct pcell *pl, struct pcell *pr)
{
    // hllc
    struct interface_val S = compute_wave_speed(pl, pr);

    struct ccell Ul = cell_primitive_to_conservative(pl);
    struct ccell Ur = cell_primitive_to_conservative(pr);
    struct fcell Fl = get_flux(pl);
    struct fcell Fr = get_flux(pr);

    real Sstar = pr->p - pl->p 
                 + Ul.ru * (S.L - pl->u) - Ur.ru * (S.R - pr->u);
    // if(Sstar != 0.0)
        Sstar /= (Ul.r * (S.L - pl->u) - Ur.r * (S.R - pr->u));
    // TODO : attention division par zéro


    struct fcell flux_out;

    if(0.0 <= S.L){
        flux_out.r  = Fl.r;
        flux_out.ru = Fl.ru;
        flux_out.e  = Fl.e;
    }
    else if(S.L <= 0.0 && 0.0 <= Sstar){
        struct fcell shift_Fl = get_flux_star_shift(&Ul, pl, S.L, Sstar);
        flux_out.r  = Fl.r  + shift_Fl.r;
        flux_out.ru = Fl.ru + shift_Fl.ru;
        flux_out.e  = Fl.e  + shift_Fl.e;
    }
    else if(Sstar <= 0.0 && 0.0 <= S.R){
        struct fcell shift_Fr = get_flux_star_shift(&Ur, pr, S.R, Sstar);
        flux_out.r  = Fr.r  + shift_Fr.r;
        flux_out.ru = Fr.ru + shift_Fr.ru;
        flux_out.e  = Fr.e  + shift_Fr.e;
    }
    else if(S.R <= 0.0){
        flux_out.r  = Fr.r;
        flux_out.ru = Fr.ru;
        flux_out.e  = Fr.e;
    }
    else{
        printf("riemman solve_fluxes if/else error\n");
        printf("S : %lf, %lf, %lf\n", S.L, Sstar, S.R);
        printf("r : %lf, %lf\n", pl->r, pr->r);
        printf("u : %lf, %lf\n", pl->u, pr->u);
        printf("p : %lf, %lf\n", pl->p, pr->p);

        real gamma = grid.gamma;

        real z  = (gamma - 1.0) / (2.0 * gamma);
        
        real aL = SQRT(pl->p * grid.gamma / pl->r);
        real aR = SQRT(pr->p * grid.gamma / pr->r);

        real numerator = (aL + aR - gamma * z * (pr->u - pl->u));
        real denominator = aL / POW(pl->p, z) + aR / POW(pr->p, z);
        real ps = POW(numerator / denominator, 1.0/z);

        printf("a : %lf, %lf\n", aL, aR);
        printf("spe : %lf, %lf\n", aL + aR, gamma * z * (pr->u - pl->u));
        printf("ps : %lf\n", ps);
        assert(false);
    }
    
    return flux_out;
}