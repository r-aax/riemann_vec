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
double g1;

/// \brief Gama 2.
double g2;

/// \brief Gama 3.
double g3;

/// \brief Gama 4.
double g4;

/// \brief Gama 5.
double g5;

/// \brief Gama 6.
double g6;

/// \brief Gama 7.
double g7;

/// \brief Gama 8.
double g8;

/// \brief Init gamas values.
void init_gamas()
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
/// TODO:
void guessp(double dl, double ul, double pl, double cl,
            double dr, double ur, double pr, double cr,
            double &pm)
{
    // purpose: to provide a guessed value for pressure
    //          pm in the Star Region. The choice is made
    //          according to adaptive Riemann solver using
    //          the PVRS, TRRS and TSRS approximate
    //          Riemann solvers. See Sect. 9.5 of Chapt. 9 of Ref. 1

    double cup, gel, ger, pmax, pmin, ppv, pq, ptl, ptr,
           qmax, quser, um;

    quser = 2.0;

    // compute guess pressure from PVRS Riemann solver
    cup = 0.25*(dl + dr)*(cl + cr);
    ppv = 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
    ppv = max(0.0, ppv);
    pmin = min(pl,  pr);
    pmax = max(pl,  pr);
    qmax = pmax/pmin;

    if (qmax <= quser && (pmin <= ppv && ppv <= pmax))
        pm = ppv;     // select PVRS Riemann solver
    else {
        if (ppv < pmin) {
            // select Two-Rarefaction Riemann solver
            pq = pow(pl/pr, g1);
            um = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
            ptl = 1.0 + g7*(ul - um)/cl;
            ptr = 1.0 + g7*(um - ur)/cr;
            pm = 0.5*(pow(pl*ptl, g3) + pow(pr*ptr, g3));
        } else {
            // select Two-Shock Riemann solver with PVRS as estimate
            gel = sqrt((g5/dl)/(g6*pl + ppv));
            ger = sqrt((g5/dr)/(g6*pr + ppv));
            pm = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
        }
    }
}

void prefun(
    double &f,
    double &fd,
    double &p,
    double &dk,
    double &pk,
    double &ck)
{
    // purpose: to evaluate the pressure functions
    //          fl and fr in exact Riemann solver
    //          and their first derivatives

    double ak, bk, pratio, qrt;

    if (p <= pk) {
        // rarefaction wave
        pratio = p/pk;
        f = g4*ck*(pow(pratio, g1) - 1.0);
        fd = (1.0/(dk*ck))*pow(pratio, -g2);
    } else {
        //  shock wave
        ak = g5/dk;
        bk = g6*pk;
        qrt = sqrt(ak/(bk + p));
        f = (p - pk)*qrt;
        fd = (1.0 - 0.5*(p - pk)/(bk + p))*qrt;
    }
}

/// \brief
///
/// TODO:
void starpu(double dl, double ul, double pl, double cl,
            double dr, double ur, double pr, double cr,
            double &p, double &u)
{
    // purpose: to compute the solution for pressure and
    //          velocity in the Star Region

    const int nriter = 20;
    const double tolpre = 1.0e-6;
    double change, fl, fld, fr, frd, pold, pstart, udiff;

    // guessed value pstart is computed
    guessp(dl, ul, pl, cl, dr, ur, pr, cr, pstart);
    pold = pstart;
    udiff = ur - ul;

    cout << "----------------------------------------\n"
         << "   Iteration number     Change\n"
         << "----------------------------------------" << endl;

    int i = 1;
    for ( ; i <= nriter; i++) {
        prefun(fl, fld, pold, dl, pl, cl);
        prefun(fr, frd, pold, dr, pr, cr);
        p = pold - (fl + fr + udiff)/(fld + frd);
        change = 2.0*abs((p - pold)/(p + pold));
        cout << '\t' << i <<  "\t\t" << change << endl;
        if (change <= tolpre)
            break;
        if (p < 0.0)
            p = tolpre;
        pold = p;
    }
    if (i > nriter)
        cout << "divergence in Newton-Raphson iteration" << endl;

    // compute velocity in star region
    u = 0.5*(ul + ur + fr - fl);
}

/// \brief
///
/// TODO:
void sample(double dl, double ul, double pl, double cl,
            double dr, double ur, double pr, double cr,
            const double pm, const double um, const double s,
            double &d, double &u, double &p)
{
    // purpose: to sample the solution throughout the wave
    //          pattern. Pressure pm and velocity um in the
    //          star region are known. Sampling is performed
    //          in terms of the 'speed' s = x/t. Sampled
    //          values are d, u, p

    double c, cml, cmr, pml, pmr, shl, shr, sl, sr, stl, str;

    if (s <= um) {
        // sampling point lies to the left of the contact discontinuity
        if (pm <= pl) {
            // left rarefaction
            shl = ul - cl;
            if (s <= shl) {
                // sampled point is left data state
                d = dl;
                u = ul;
                p = pl;
            } else {
                cml = cl*pow(pm/pl, g1);
                stl = um - cml;
                if (s > stl) {
                    // sampled point is star left state
                    d = dl*pow(pm/pl, 1.0/GAMA);
                    u = um;
                    p = pm;
                } else {
                    // sampled point is inside left fan
                    u = g5*(cl + g7*ul + s);
                    c = g5*(cl + g7*(ul - s));
                    d = dl*pow(c/cl, g4);
                    p = pl*pow(c/cl, g3);
                }
            }
        } else {
            // left shock
            pml = pm/pl;
            sl = ul - cl*sqrt(g2*pml + g1);
            if (s <= sl) {
                // sampled point is left data state
                d = dl;
                u = ul;
                p = pl;
            } else {
                // sampled point is star left state
                d = dl*(pml + g6)/(pml*g6 + 1.0);
                u = um;
                p = pm;
            }
        }
    } else {
        // sampling point lies to the right of the contact discontinuity
        if (pm > pr) {
            // right shock
            pmr = pm/pr;
            sr  = ur + cr*sqrt(g2*pmr + g1);
            if (s >= sr) {
                // sampled point is right data state
                d = dr;
                u = ur;
                p = pr;
            } else {
                // sampled point is star right state
                d = dr*(pmr + g6)/(pmr*g6 + 1.0);
                u = um;
                p = pm;
            }
        } else {
            // right rarefaction
            shr = ur + cr;
            if (s >= shr) {
                // sampled point is right data state
                d = dr;
                u = ur;
                p = pr;
            } else {
                cmr = cr*pow(pm/pr, g1);
                str = um + cmr;
                if (s <= str) {
                    // sampled point is star right state
                    d = dr*pow(pm/pr, 1.0/GAMA);
                    u = um;
                    p = pm;
                } else {
                    // sampled point is inside left fan
                    u = g5*(-cr + g7*ur + s);
                    c = g5*(cr - g7*(ur - s));
                    d = dr*pow(c/cr, g4);
                    p = pr*pow(c/cr, g3);
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
void riemann(double dl, double ul, double pl,
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

int main()
{
    return 0;
}
