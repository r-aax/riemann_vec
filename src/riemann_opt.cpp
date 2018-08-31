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
#include "riemann.h"
using namespace std;

#ifdef INTEL

#include <immintrin.h>

/// \brief Short name for LD intrinsic.
#define LD(ADDR) _mm512_load_ps(ADDR)

/// \brief Short name for ST intrinsic.
#define ST(ADDR, VAL) _mm512_store_ps(ADDR, VAL)

/// \brief Short name for ADD intrinsic.
#define ADD(va, vb) _mm512_add_ps(va, vb)

/// \brief Short name for MUL intrinsic.
#define MUL(va, vb) _mm512_mul_ps(va, vb)

/// \brief Short name for SUB intrinsic.
#define SUB(va, vb) _mm512_sub_ps(va, vb)

/// \brief Short name for DIV intrinsic.
#define DIV(va, vb) _mm512_div_ps(va, vb)

/// \brief Short name for POW intrinsic.
#define POW(va, vb) _mm512_pow_ps(va, vb)

/// \brief Short name for SQRT intrinsic.
#define SQRT(va) _mm512_sqrt_ps(va)

/// \brief Zero.
__m512 z = _mm512_setzero_ps();

/// \brief 1.
__m512 v1 = _mm512_set1_ps(1.0);

/// \brief G1.
__m512 g1 = _mm512_set1_ps(G1);

/// \brief G2.
__m512 g2 = _mm512_set1_ps(G2);

/// \brief G3.
__m512 g3 = _mm512_set1_ps(G3);

/// \brief G4.
__m512 g4 = _mm512_set1_ps(G4);

/// \brief G5.
__m512 g5 = _mm512_set1_ps(G5);

/// \brief G6.
__m512 g6 = _mm512_set1_ps(G6);

/// \brief G7.
__m512 g7 = _mm512_set1_ps(G7);

/// \brief 1/GAMA.
__m512 igama = _mm512_set1_ps(1.0 / GAMA);

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
static void guessp(float dl, float ul, float pl, float cl,
                   float dr, float ur, float pr, float cr,
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

/// \brief
///
/// Purpose is to evaluate the pressure functions
/// fl and fr in exact Riemann solver
/// and their first derivatives.
static void prefun(float &f, float &fd, float &p,
                   float &dk, float &pk, float &ck)
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

/// \brief
///
/// Purpose is to compute the solution for pressure
/// and velocity in the Star Region.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound velocity
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound velocity
/// \param[in] p - pressure in star region
/// \param[out] u - velocity in star region
static void starpu(float dl, float ul, float pl, float cl,
                   float dr, float ur, float pr, float cr,
                   float &p, float &u)
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
    u = 0.5*(ul + ur + fr - fl);
}

/// \brief
///
/// Purpose is to compute the solution for pressure
/// and velocity in the Star Region.
///
/// \param[in] vdl - left side density
/// \param[in] vul - left side velocity
/// \param[in] vpl - left side pressure
/// \param[in] vcl - left side sound velocity
/// \param[in] vdr - right side density
/// \param[in] vur - right side velocity
/// \param[in] vpr - right side pressure
/// \param[in] vcr - right side sound velocity
/// \param[out] vp - pressure in star region
/// \param[out] vu - velocity in star region
static void starpu_16(__m512 vdl, __m512 vul, __m512 vpl, __m512 vcl,
                      __m512 vdr, __m512 vur, __m512 vpr, __m512 vcr,
                      __m512 *vp, __m512 *vu)
{
    float arr_dl[16], arr_ul[16], arr_pl[16], arr_cl[16],
          arr_dr[16], arr_ur[16], arr_pr[16], arr_cr[16],
          arr_p[16], arr_u[16];

    _mm512_store_ps(&arr_dl[0], vdl);
    _mm512_store_ps(&arr_ul[0], vul);
    _mm512_store_ps(&arr_pl[0], vpl);
    _mm512_store_ps(&arr_cl[0], vcl);
    _mm512_store_ps(&arr_dr[0], vdr);
    _mm512_store_ps(&arr_ur[0], vur);
    _mm512_store_ps(&arr_pr[0], vpr);
    _mm512_store_ps(&arr_cr[0], vcr);

    for (int i = 0; i < 16; i++)
    {
        const int nriter = 20;
        const float tolpre = 1.0e-6;
        float change, fl, fld, fr, frd, pold, pstart, udiff;

        // Guessed value pstart is computed.
        guessp(arr_dl[i], arr_ul[i], arr_pl[i], arr_cl[i], arr_dr[i], arr_ur[i], arr_pr[i], arr_cr[i], pstart);
        pold = pstart;
        udiff = arr_ur[i] - arr_ul[i];

        int i = 1;

        for ( ; i <= nriter; i++)
        {
            prefun(fl, fld, pold, arr_dl[i], arr_pl[i], arr_cl[i]);
            prefun(fr, frd, pold, arr_dr[i], arr_pr[i], arr_cr[i]);
            arr_p[i] = pold - (fl + fr + udiff) / (fld + frd);
            change = 2.0 * abs((arr_p[i] - pold) / (arr_p[i] + pold));

            if (change <= tolpre)
            {
                break;
            }

            if (arr_p[i] < 0.0)
            {
                arr_p[i] = tolpre;
            }

            pold = arr_p[i];
        }

        if (i > nriter)
        {
            cout << "divergence in Newton-Raphson iteration" << endl;

            exit(1);
        }

        // compute velocity in star region
        arr_u[i] = 0.5*(arr_ul[i] + arr_ur[i] + fr - fl);
    }

    *vp = _mm512_load_ps(&arr_p[0]);
    *vu = _mm512_load_ps(&arr_u[0]);
}

/// \brief
///
/// Purpose is to sample the solution throughout the wave
/// pattern. Pressure pm and velocit
/// star region are known. Sampling is performed
/// in terms of the 'speed' s = x/t. Sampled
/// values are d, u, p.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound velocity
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound velocity
/// \param[out] od - result density
/// \param[out] ou - result velocity
/// \param[out] op - result pressure
static void sample_16(__m512 dl, __m512 ul, __m512 pl, __m512 cl,
                      __m512 dr, __m512 ur, __m512 pr, __m512 cr,
                      __m512 pm, __m512 um,
                      __m512 *od, __m512 *ou, __m512 *op)
{
    __m512 d, u, p, c, ums, pms, sh, st, s, uc;
    __mmask16 cond_um, cond_pm, cond_sh, cond_st, cond_s, cond_sh_st;

    // d/u/p/c/ums
    cond_um = _mm512_cmp_ps_mask(um, z, _MM_CMPINT_LT);
    d = _mm512_mask_blend_ps(cond_um, dl, dr);
    u = _mm512_mask_blend_ps(cond_um, ul, ur);
    p = _mm512_mask_blend_ps(cond_um, pl, pr);
    c = _mm512_mask_blend_ps(cond_um, cl, cr);
    ums = um;
    u = _mm512_mask_sub_ps(u, cond_um, z, u);
    ums = _mm512_mask_sub_ps(ums, cond_um, z, ums);

    // Calculate main values.
    pms = DIV(pm, p);
    sh = SUB(u, c);
    st = _mm512_fnmadd_ps(POW(pms, g1), c, ums);
    s = _mm512_fnmadd_ps(c, SQRT(_mm512_fmadd_ps(g2, pms, g1)), u);

    // Conditions.
    cond_pm = _mm512_cmp_ps_mask(pm, p, _MM_CMPINT_LE);
    cond_sh = _mm512_mask_cmp_ps_mask(cond_pm, sh, z, _MM_CMPINT_LT);
    cond_st = _mm512_mask_cmp_ps_mask(cond_sh, st, z, _MM_CMPINT_LT);
    cond_s = _mm512_mask_cmp_ps_mask(~cond_pm, s, z, _MM_CMPINT_LT);

    // Store.
    d = _mm512_mask_mov_ps(d, cond_st, MUL(d, POW(pms, igama)));
    d = _mm512_mask_mov_ps(d, cond_s, MUL(d, DIV(ADD(pms, g6), _mm512_fmadd_ps(pms, g6, v1))));
    u = _mm512_mask_mov_ps(u, cond_st | cond_s, ums);
    p = _mm512_mask_mov_ps(p, cond_st | cond_s, pm);

    // Low prob - ignnore it.
    cond_sh_st = cond_sh & ~cond_st;
    if (cond_sh_st)
    {
        u = _mm512_mask_mov_ps(u, cond_sh_st, MUL(g5, _mm512_fmadd_ps(g7, u, c)));
        uc = DIV(u, c);
        d = _mm512_mask_mov_ps(d, cond_sh_st, MUL(d, POW(uc, g4)));
        p = _mm512_mask_mov_ps(p, cond_sh_st, MUL(p, POW(uc, g3)));
    }

    // Final store.
    u = _mm512_mask_sub_ps(u, cond_um, z, u);
    *od = d;
    *ou = u;
    *op = p;
}

/// \brief Riemann solver for 16 cases.
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
static void riemann_16(float *dl, float *ul, float *pl,
                       float *dr, float *ur, float *pr,
                       float *d, float *u, float *p)
{
    __assume_aligned(dl, 64);
    __assume_aligned(ul, 64);
    __assume_aligned(pl, 64);
    __assume_aligned(dr, 64);
    __assume_aligned(ur, 64);
    __assume_aligned(pr, 64);
    __assume_aligned(d, 64);
    __assume_aligned(u, 64);
    __assume_aligned(p, 64);

    // Basic data.
    __m512 vdl = _mm512_load_ps(dl);
    __m512 vul = _mm512_load_ps(ul);
    __m512 vpl = _mm512_load_ps(pl);
    __m512 vdr = _mm512_load_ps(dr);
    __m512 vur = _mm512_load_ps(ur);
    __m512 vpr = _mm512_load_ps(pr);

    // Sound speed and check for vacuum.
    __m512 vgama = _mm512_set1_ps(GAMA);
    __m512 vg4 = _mm512_set1_ps(G4);
    __m512 vcl = _mm512_sqrt_ps(_mm512_div_ps(_mm512_mul_ps(vgama, vpl), vdl));
    __m512 vcr = _mm512_sqrt_ps(_mm512_div_ps(_mm512_mul_ps(vgama, vpr), vdr));
    __mmask16 vacuum_mask = _mm512_cmp_ps_mask(_mm512_mul_ps(vg4, _mm512_add_ps(vcl, vcr)),
                                               _mm512_sub_ps(vur, vul),
                                               _MM_CMPINT_LE);

    if (vacuum_mask != 0x0)
    {
        cerr << "VACUUM" << endl;
        exit(1);
    }

    __m512 vpm, vum, vd, vu, vp;
    starpu_16(vdl, vul, vpl, vcl, vdr, vur, vpr, vcr, &vpm, &vum);
    sample_16(vdl, vul, vpl, vcl, vdr, vur, vpr, vcr, vpm, vum, &vd, &vu, &vp);
    _mm512_store_ps(d, vd);
    _mm512_store_ps(u, vu);
    _mm512_store_ps(p, vp);
}

#endif

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
                 float *dl, float *ul, float *pl,
                 float *dr, float *ur, float *pr,
                 float *d, float *u, float *p)
{

#ifndef INTEL

    riemann(c, dl, ul, pl, dr, ur, pr, d, u, p);

#else

    float d_, u_, p_;

    int c_tail = c & 0xF;
    int c_without_tail = c - c_tail;

    // Main body.
    for (int i = 0; i < c_without_tail; i += 16)
    {
        riemann_16(&dl[i], &ul[i], &pl[i], &dr[i], &ur[i], &pr[i], &d[i], &u[i], &p[i]);
    }

    // Tail.
    for (int i = c_without_tail; i < c; i++)
    {
        riemann(dl[i], ul[i], pl[i], dr[i], ur[i], pr[i], d_, u_, p_);
        d[i] = d_;
        u[i] = u_;
        p[i] = p_;
    }

#endif

}
