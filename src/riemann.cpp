/// \file
/// \brief Riemann solver.
///
/// Exact Riemann solver for the Euler equations in one dimension
/// Translated from the Fortran code er1pex.f and er1pex.ini
/// by Dr. E.F. Toro downloaded from
/// http://www.numeritek.com/numerica_software.html#freesample

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <omp.h>
#include "riemann.h"

using namespace std;

/// \brief Calculate start guessed pressure value.
///
/// Purpose is to provide a guessed value for pressure
/// pm in the Star Region. The choice is made
/// according to adaptive Riemann solver using
/// the PVRS, TRRS and TSRS approximate
/// Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound speed
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound speed
/// \param[out] pm - pressure
static void
guessp(float dl,
       float ul,
       float pl,
       float cl,
       float dr,
       float ur,
       float pr,
       float cr,
       float &pm)
{
    float cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr, qmax, quser, um;

    quser = 2.0;

    // Compute guess pressure from PVRS Riemann solver.
    cup = 0.25 * (dl + dr) * (cl + cr);
    ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
    ppv = (ppv > 0.0) ? ppv : 0.0;
    pmin = (pl < pr) ? pl : pr;
    pmax = (pl > pr) ? pl : pr;
    qmax = pmax / pmin;

    if ((qmax <= quser) && (pmin <= ppv) && (ppv <= pmax))
    {
        // Select PVRS Riemann solver.
        pm = ppv;
    }
    else
    {
        if (ppv < pmin)
        {
            // Select Two-Rarefaction Riemann solver.
            pq = pow(pl / pr, G1);
            um = (pq * ul / cl + ur / cr + G4 * (pq - 1.0)) / (pq / cl + 1.0 / cr);
            ptl = 1.0 + G7 * (ul - um) / cl;
            ptr = 1.0 + G7 * (um - ur) / cr;
            pm = 0.5 * (pow(pl * ptl, G3) + pow(pr * ptr, G3));
        }
        else
        {
            // Select Two-Shock Riemann solver with PVRS as estimate.
            gel = sqrt((G5 / dl) / (G6 * pl + ppv));
            ger = sqrt((G5 / dr) / (G6 * pr + ppv));
            pm = (gel * pl + ger * pr - (ur - ul)) / (gel + ger);
        }
    }
}

/// \brief Calculate pressures and derivatives for left and right sides.
///
/// Purpose is to evaluate the pressure functions
/// fl and fr in exact Riemann solver
/// and their first derivatives.
///
/// \param[out] f - pressure function
/// \param[out] fd - pressure derivative
/// \param[in] p - old pressure
/// \param[in] dk - density
/// \param[in] pk - pressure
/// \param[in] ck - sound speed
static void
prefun(float &f,
       float &fd,
       float &p,
       float &dk,
       float &pk,
       float &ck)
{
    float ak, bk, pratio, qrt;

    if (p <= pk)
    {
        // Rarefaction wave.
        pratio = p / pk;
        f = G4 * ck * (pow(pratio, G1) - 1.0);
        fd = (1.0 / (dk * ck)) * pow(pratio, -G2);
    }
    else
    {
        // Shock wave.
        ak = G5 / dk;
        bk = G6 * pk;
        qrt = sqrt(ak / (bk + p));
        f = (p - pk) * qrt;
        fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
    }
}

/// \brief Pressure and speed calculation in star region.
///
/// Purpose is to compute the solution for pressure
/// and velocity in the Star Region.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side speed
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound speed
/// \param[in] dr - right side density
/// \param[in] ur - right side speed
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound speed
/// \param[out] p - pressure in star region
/// \param[out] u - speed in star region
static void
starpu(float dl,
       float ul,
       float pl,
       float cl,
       float dr,
       float ur,
       float pr,
       float cr,
       float &p,
       float &u)
{
    const int nriter = 20;
    const float tolpre = 1.0e-6;
    float change, fl, fld, fr, frd, pold, pstart, udiff;

    // Guessed value pstart is computed.
    guessp(dl, ul, pl, cl, dr, ur, pr, cr, pstart);
    pold = pstart;
    udiff = ur - ul;

    int i = 1;

    for ( ; i <= nriter; i++)
    {
        prefun(fl, fld, pold, dl, pl, cl);
        prefun(fr, frd, pold, dr, pr, cr);
        p = pold - (fl + fr + udiff) / (fld + frd);
        change = 2.0 * abs((p - pold) / (p + pold));

        if (change <= tolpre)
        {
            break;
        }

        if (p < 0.0)
        {
            p = tolpre;
        }

        pold = p;
    }

    if (i > nriter)
    {
        cout << "divergence in Newton-Raphson iteration" << endl;

        exit(1);
    }

    // compute velocity in star region
    u = 0.5 * (ul + ur + fr - fl);
}

/// \brief Final analyze of the configuration.
///
/// Purpose is to sample the solution throughout the wave
/// pattern. Pressure pm and velocity um in the
/// star region are known. Sampling is performed
/// in terms of the 'speed' s = x/t. Sampled
/// values are d, u, p.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side speed x component
/// \param[in] vl - left side speed y component
/// \param[in] wl - left side speed z component
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound speed
/// \param[in] dr - right side density
/// \param[in] ur - right side speed x component
/// \param[in] vr - right side speed y component
/// \param[in] wr - right side speed z component
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound speed
/// \param[in] pm - pressure in star region
/// \param[in] um - speed in star region
/// \param[out] d - density
/// \param[out] u - speed x component
/// \param[out] v - speed y component
/// \param[out] w - speed z component
/// \param[out] p - pressure
static void
sample(float dl,
       float ul,
       float vl,
       float wl,
       float pl,
       float cl,
       float dr,
       float ur,
       float vr,
       float wr,
       float pr,
       float cr,
       const float pm,
       const float um,
       float &d,
       float &u,
       float &v,
       float &w,
       float &p)
{
    float c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

    if (0.0 <= um)
    {
        // Sampling point lies to the left of the contact discontinuity.
        v = vl;
        w = wl;

        if (pm <= pl)
        {
            // Left rarefaction.
            shl = ul - cl;

            if (0.0 <= shl)
            {
                // Sampled point is left data state.
                d = dl;
                u = ul;
                p = pl;
            }
            else
            {
                cml = cl * pow(pm / pl, G1);
                stl = um - cml;

                if (0.0 > stl)
                {
                    // Sampled point is star left state.
                    d = dl * pow(pm / pl, 1.0 / GAMA);
                    u = um;
                    p = pm;
                }
                else
                {
                    // Sampled point is inside left fan.
                    u = G5 * (cl + G7 * ul);
                    c = G5 * (cl + G7 * ul);
                    d = dl * pow(c / cl, G4);
                    p = pl * pow(c / cl, G3);
                }
            }
        }
        else
        {
            // Left shock.
            pml = pm / pl;
            sl = ul - cl * sqrt(G2 * pml + G1);

            if (0.0 <= sl)
            {
                // Sampled point is left data state.
                d = dl;
                u = ul;
                p = pl;
            }
            else
            {
                // Sampled point is star left state.
                d = dl * (pml + G6) / (pml * G6 + 1.0);
                u = um;
                p = pm;
            }
        }
    }
    else
    {
        // Sampling point lies to the right of the contact discontinuity.
        v = vr;
        w = wr;

        if (pm > pr)
        {
            // Right shock.
            pmr = pm / pr;
            sr  = ur + cr * sqrt(G2 * pmr + G1);

            if (0.0 >= sr)
            {
                // Sampled point is right data state.
                d = dr;
                u = ur;
                p = pr;
            }
            else
            {
                // Sampled point is star right state.
                d = dr * (pmr + G6) / (pmr * G6 + 1.0);
                u = um;
                p = pm;
            }
        }
        else
        {
            // Right rarefaction.
            shr = ur + cr;
            if (0.0 >= shr)
            {
                // Sampled point is right data state.
                d = dr;
                u = ur;
                p = pr;
            }
            else
            {
                cmr = cr * pow(pm / pr, G1);
                str = um + cmr;

                if (0.0 <= str)
                {
                    // Sampled point is star right state.
                    d = dr * pow(pm / pr, 1.0 / GAMA);
                    u = um;
                    p = pm;
                }
                else
                {
                    // Sampled point is inside left fan.
                    u = G5 * (-cr + G7 * ur);
                    c = G5 * (cr - G7 * ur);
                    d = dr * pow(c / cr, G4);
                    p = pr * pow(c / cr, G3);
                }
            }
        }
    }
}

/// \brief Riemann solver.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference x component
/// \param[out] v - result velocity reference y component
/// \param[out] w - result velocity reference z component
/// \param[out] p - result pressure reference
static void
riemann(float dl,
        float ul,
        float vl,
        float wl,
        float pl,
        float dr,
        float ur,
        float vr,
        float wr,
        float pr,
        float &d,
        float &u,
        float &v,
        float &w,
        float &p)
{
    float pm, um, cl, cr;

    pm = 0.0;

    // Sound speeds.
    cl = sqrt(GAMA * pl / dl);
    cr = sqrt(GAMA * pr / dr);

    // Check for vacuum.
    if (G4 * (cl + cr) <= (ur - ul))
    {
        cerr << "VACUUM" << endl;
        exit(1);
    }

    // Exact solution.
    starpu(dl, ul, pl, cl, dr, ur, pr, cr, pm, um);
    sample(dl, ul, vl, wl, pl, cl,
           dr, ur, vr, wr, pr, cr,
           pm, um,
           d, u, v, w, p);
}

/// \brief Riemann not optimize for multiple data.
///
/// \param[in] c - sizes
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference x component
/// \param[out] v - result velocity reference y component
/// \param[out] w - result velocity reference z component
/// \param[out] p - result pressure reference
void
riemann_s(int c,
          float *dl,
          float *ul,
          float *vl,
          float *wl,
          float *pl,
          float *dr,
          float *ur,
          float *vr,
          float *wr,
          float *pr,
          float *d,
          float *u,
          float *v,
          float *w,
          float *p)
{
    float d_, u_, v_, w_, p_;

    for (int i = 0; i < c; i++)
    {
        riemann(dl[i], ul[i], vl[i], wl[i], pl[i],
                dr[i], ur[i], vr[i], wr[i], pr[i],
                d_, u_, v_, w_, p_);
        d[i] = d_;
        u[i] = u_;
        v[i] = v_;
        w[i] = w_;
        p[i] = p_;
    }
}

/// \brief Riemann 16x not optimized solver.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference x component
/// \param[out] v - result velocity reference y component
/// \param[out] w - result velocity reference z component
/// \param[out] p - result pressure reference
void
riemann_16_s(float *dl,
             float *ul,
             float *vl,
             float *wl,
             float *pl,
             float *dr,
             float *ur,
             float *vr,
             float *wr,
             float *pr,
             float *d,
             float *u,
             float *v,
             float *w,
             float *p)
{
    riemann_s(16,
              dl, ul, vl, wl, pl,
              dr, ur, vr, wr, pr,
              d, u, v, w, p);
}

/// \brief Riemann solver for multiple data.
///
/// \param[in] c - cases count
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference x component
/// \param[out] v - result velocity reference y component
/// \param[out] w - result velocity reference z component
/// \param[out] p - result pressure reference
/// \param[in] - solver function fo 16x data
void
riemann_n(int c,
          float *dl,
          float *ul,
          float *vl,
          float *wl,
          float *pl,
          float *dr,
          float *ur,
          float *vr,
          float *wr,
          float *pr,
          float *d,
          float *u,
          float *v,
          float *w,
          float *p,
          int nt,
          void (*solver_16)(float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *,
                            float *))
{
    int c_tail = c & 0xF;
    int c_base = c - c_tail;

    //
    // Main body.
    //

    omp_set_num_threads(nt);

    #pragma omp parallel
    {
        int tn = omp_get_thread_num();
        int lb = (int)((c / FP16_VECTOR_SIZE) * ((double)tn / (double)nt));
        int ub = (int)((c / FP16_VECTOR_SIZE) * ((double)(tn + 1) / (double)nt));

        for (int i = lb * FP16_VECTOR_SIZE;
             i < ub * FP16_VECTOR_SIZE;
             i += FP16_VECTOR_SIZE)
        {
            solver_16(dl + i, ul + i, vl + i, wl + i, pl + i,
                      dr + i, ur + i, vr + i, wr + i, pr + i,
                      d + i, u + i, v + i, w + i, p + i);
        }
    }

    omp_set_num_threads(1);

    // Process tail.
    riemann_s(c_tail,
              dl + c_base, ul + c_base, vl + c_base, wl + c_base, pl + c_base,
              dr + c_base, ur + c_base, vr + c_base, wr + c_base, pr + c_base,
              d + c_base, u + c_base, v + c_base, w + c_base, p + c_base);
}

/// \brief Riemann solver not vectorized for scalar data.
///
/// \param[in] c - cases count
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference x component
/// \param[out] v - result velocity reference y component
/// \param[out] w - result velocity reference z component
/// \param[out] p - result pressure reference
/// \param[out] nt - number of threads
void
riemann_n_s(int c,
            float *dl,
            float *ul,
            float *vl,
            float *wl,
            float *pl,
            float *dr,
            float *ur,
            float *vr,
            float *wr,
            float *pr,
            float *d,
            float *u,
            float *v,
            float *w,
            float *p,
            int nt)
{
    riemann_n(c,
              dl, ul, vl, wl, pl,
              dr, ur, vr, wr, pr,
              d, u, v, w, p,
              nt,
              riemann_16_s);
}
