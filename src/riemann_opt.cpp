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
using namespace std;

/// \brief Gama value.
#define GAMA 1.4

/// \brief Gama 1.
static double g1;

/// \brief Gama 2.
static double g2;

/// \brief Gama 3.
static double g3;

/// \brief Gama 4.
static double g4;

/// \brief Gama 5.
static double g5;

/// \brief Gama 6.
static double g6;

/// \brief Gama 7.
static double g7;

/// \brief Gama 8.
static double g8;

/// \brief Init gamas values.
void init_gamas_opt()
{
    g1 = (GAMA - 1.0) / (2.0 * GAMA);
    g2 = (GAMA + 1.0) / (2.0 * GAMA);
    g3 = 2.0 * GAMA / (GAMA - 1.0);
    g4 = 2.0 / (GAMA - 1.0);
    g5 = 2.0 / (GAMA + 1.0);
    g6 = (GAMA - 1.0) / (GAMA + 1.0);
    g7 = (GAMA - 1.0) / 2.0;
    g8 = GAMA - 1.0;
}

/// \brief 
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
static void guessp(double dl, double ul, double pl, double cl,
                   double dr, double ur, double pr, double cr,
                   double &pm)
{
    double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr, qmax, quser, um;

    quser = 2.0;

    // Compute guess pressure from PVRS Riemann solver.
    cup = 0.25 * (dl + dr) * (cl + cr);
    ppv = 0.5 * (pl + pr) + 0.5 * (ul - ur) * cup;
    ppv = max(0.0, ppv);
    pmin = min(pl,  pr);
    pmax = max(pl,  pr);
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
            pq = pow(pl / pr, g1);
            um = (pq * ul / cl + ur / cr + g4 * (pq - 1.0)) / (pq / cl + 1.0 / cr);
            ptl = 1.0 + g7 * (ul - um) / cl;
            ptr = 1.0 + g7 * (um - ur) / cr;
            pm = 0.5 * (pow(pl * ptl, g3) + pow(pr * ptr, g3));
        }
        else
        {
            // Select Two-Shock Riemann solver with PVRS as estimate.
            gel = sqrt((g5 / dl) / (g6 * pl + ppv));
            ger = sqrt((g5 / dr) / (g6 * pr + ppv));
            pm = (gel * pl + ger * pr - (ur - ul)) / (gel + ger);
        }
    }
}

/// \brief
///
/// Purpose is to evaluate the pressure functions
/// fl and fr in exact Riemann solver
/// and their first derivatives.
///
/// TODO:
static void prefun(double &f, double &fd, double &p,
                   double &dk, double &pk, double &ck)
{
    double ak, bk, pratio, qrt;

    if (p <= pk)
    {
        // Rarefaction wave.
        pratio = p / pk;
        f = g4 * ck * (pow(pratio, g1) - 1.0);
        fd = (1.0 / (dk * ck)) * pow(pratio, -g2);
    }
    else
    {
        // Shock wave.
        ak = g5 / dk;
        bk = g6 * pk;
        qrt = sqrt(ak / (bk + p));
        f = (p - pk) * qrt;
        fd = (1.0 - 0.5 * (p - pk) / (bk + p)) * qrt;
    }
}

/// \brief
///
/// Purpose is to compute the solution for pressure
/// and velocity in the Star Region.
///
/// TODO:
static void starpu(double dl, double ul, double pl, double cl,
                   double dr, double ur, double pr, double cr,
                   double &p, double &u)
{
    const int nriter = 20;
    const double tolpre = 1.0e-6;
    double change, fl, fld, fr, frd, pold, pstart, udiff;

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
    u = 0.5*(ul + ur + fr - fl);
}

/// \brief
///
/// Purpose is to sample the solution throughout the wave
/// pattern. Pressure pm and velocity um in the
/// star region are known. Sampling is performed
/// in terms of the 'speed' s = x/t. Sampled
/// values are d, u, p.
///
/// TODO:
static void sample(double dl, double ul, double pl, double cl,
                   double dr, double ur, double pr, double cr,
                   const double pm, const double um, const double s,
                   double &d, double &u, double &p)
{
    double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

    if (s <= um)
    {
        // Sampling point lies to the left of the contact discontinuity.
        if (pm <= pl)
        {
            // Left rarefaction.
            shl = ul - cl;

            if (s <= shl)
            {
                // Sampled point is left data state.
                d = dl;
                u = ul;
                p = pl;
            }
            else
            {
                cml = cl * pow(pm / pl, g1);
                stl = um - cml;

                if (s > stl)
                {
                    // Sampled point is star left state.
                    d = dl * pow(pm / pl, 1.0 / GAMA);
                    u = um;
                    p = pm;
                }
                else
                {
                    // Sampled point is inside left fan.
                    u = g5 * (cl + g7 * ul + s);
                    c = g5 * (cl + g7 * (ul - s));
                    d = dl * pow(c / cl, g4);
                    p = pl * pow(c / cl, g3);
                }
            }
        }
        else
        {
            // Left shock.
            pml = pm / pl;
            sl = ul - cl * sqrt(g2 * pml + g1);

            if (s <= sl)
            {
                // Sampled point is left data state.
                d = dl;
                u = ul;
                p = pl;
            }
            else
            {
                // Sampled point is star left state.
                d = dl * (pml + g6) / (pml * g6 + 1.0);
                u = um;
                p = pm;
            }
        }
    }
    else
    {
        // Sampling point lies to the right of the contact discontinuity.
        if (pm > pr)
        {
            // Right shock.
            pmr = pm / pr;
            sr  = ur + cr * sqrt(g2 * pmr + g1);

            if (s >= sr)
            {
                // Sampled point is right data state.
                d = dr;
                u = ur;
                p = pr;
            }
            else
            {
                // Sampled point is star right state.
                d = dr * (pmr + g6) / (pmr * g6 + 1.0);
                u = um;
                p = pm;
            }
        }
        else
        {
            // Right rarefaction.
            shr = ur + cr;
            if (s >= shr)
            {
                // Sampled point is right data state.
                d = dr;
                u = ur;
                p = pr;
            }
            else
            {
                cmr = cr * pow(pm / pr, g1);
                str = um + cmr;

                if (s <= str)
                {
                    // Sampled point is star right state.
                    d = dr * pow(pm / pr, 1.0 / GAMA);
                    u = um;
                    p = pm;
                }
                else
                {
                    // Sampled point is inside left fan.
                    u = g5 * (-cr + g7 * ur + s);
                    c = g5 * (cr - g7 * (ur - s));
                    d = dr * pow(c / cr, g4);
                    p = pr * pow(c / cr, g3);
                }
            }
        }
    }
}

/// \brief Riemann solver.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference
/// \param[out] p - result pressure reference
static void riemann(double dl, double ul, double pl,
                    double dr, double ur, double pr,
                    double &d, double &u, double &p)
{
    double pm, um, cl, cr;

    pm = 0.0;

    // Sound speeds.
    cl = sqrt(GAMA * pl / dl);
    cr = sqrt(GAMA * pr / dr);

    // Check for vacuum.
    if (g4 * (cl + cr) <= (ur - ul))
    {

        cerr << "VACUUM" << endl;
        exit(1);
    }    

    // Exact solution.
    starpu(dl, ul, pl, cl, dr, ur, pr, cr, pm, um);
    sample(dl, ul, pl, cl, dr, ur, pr, cr, pm, um, 0.0, d, u, p);
}

/// \brief Riemann solver.
///
/// \param[in] c - cases count
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference
/// \param[out] p - result pressure reference
void riemann_opt(int c,
                 double *dl, double *ul, double *pl,
                 double *dr, double *ur, double *pr,
                 double *d, double *u, double *p)
{
    double d_, u_, p_;

    for (int i = 0; i < c; i++)
    {
        riemann(dl[i], ul[i], pl[i], dr[i], ur[i], pr[i], d_, u_, p_);
        d[i] = d_;
        u[i] = u_;
        p[i] = p_;
    }
}
