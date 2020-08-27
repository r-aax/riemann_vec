/// file
/// \brief Riemann solver (vectorized version).
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

/// \brief Fp32 vector length.
#define FP16_VECTOR_SIZE 16

/// \brief Short name for load intrinsic.
#define LD(ADDR) _mm512_load_ps(ADDR)

/// \brief Short name for store intrinsic.
#define ST(ADDR, VAL) _mm512_store_ps(ADDR, VAL)

/// \brief Short name for add intrinsic.
#define ADD(VA, VB) _mm512_add_ps(VA, VB)

/// \brief Short name for mul intrinsic.
#define MUL(VA, VB) _mm512_mul_ps(VA, VB)

/// \brief Short name for sub intrinsic.
#define SUB(VA, VB) _mm512_sub_ps(VA, VB)

/// \brief Short name for div intrinsic.
#define DIV(VA, VB) _mm512_div_ps(VA, VB)

/// \brief Short name for pow intrinsic.
#define POW(VA, VB) _mm512_pow_ps(VA, VB)

/// \brief Short name for sqrt intrinsic.
#define SQRT(V) _mm512_sqrt_ps(V)

/// \brief Short name for abs intrinsic.
#define ABS(V) _mm512_abs_ps(V)

/// \brief Short name for max intrinsic.
#define MAX(VA, VB) _mm512_max_ps(VA, VB)

/// \brief Short name for min intrinsic.
#define MIN(VA, VB) _mm512_min_ps(VA, VB)

/// \bries Short name for cmp.
#define CMP(VA, VB, CMP) _mm512_cmp_ps_mask(VA, VB, CMP)

/// \brief Short name for setzero intrinsic.
#define SETZERO() _mm512_setzero_ps()

/// \brief Short name for set1 intrinsic.
#define SET1(V) _mm512_set1_ps(V)

/// \brief Short name for fmadd intrinsic.
#define FMADD(VA, VB, VC) _mm512_fmadd_ps(VA, VB, VC)

/// \brief Zero.
__m512 z = SETZERO();

/// \brief 1.
__m512 one = SET1(1.0);

/// \brief GAMA.
__m512 gama = SET1(GAMA);

/// \brief Gamma special value 1 vector.
__m512 g1 = SET1(G1);

/// \brief Gamma special value 2 vector.
__m512 g2 = SET1(G2);

/// \brief Gamma special value 3 vector.
__m512 g3 = SET1(G3);

/// \brief Gamma special value 4 vector.
__m512 g4 = SET1(G4);

/// \brief Gamma special value 5 vector.
__m512 g5 = SET1(G5);

/// \brief Gamma special value 6 vector.
__m512 g6 = SET1(G6);

/// \brief Gamma special value 7 vector.
__m512 g7 = SET1(G7);

/// \brief Get <c>i</c>-th element from vector.
///
/// \param[in] v - vector
/// \param[in] i - index
///
/// \return
/// Element.
static float
Get(__m512 v,
    int i)
{
    float arr[FP16_VECTOR_SIZE];

    ST(&arr[0], v);

    return arr[i];
}

/// \brief Set <c>i</c>-th element in vector.
///
/// \param[in,out] v - vector
/// \param[in] i - index
/// \param[in] f - value
static void
Set(__m512 *v,
    int i,
    float f)
{
    float arr[FP16_VECTOR_SIZE];

    ST(&arr[0], *v);
    arr[i] = f;
    *v = LD(&arr[0]);
}

/// \brief Print.
///
/// \param[in] v - vector
static void
Print(__m512 v)
{
    float arr[FP16_VECTOR_SIZE];

    ST(&arr[0], v);

    printf("[");

    for (int i = 0; i < FP16_VECTOR_SIZE; i++)
    {
        printf(" %f ", arr[i]);
    }

    printf("]\n");
}

/// \brief Count of '1' bits in mask.
///
/// \param[in] mask - mask
///
/// \return
/// Count of '1' bits in mask.
static int
cnt(int mask)
{
    int c = 0;

    for (int i = 0; i < FP16_VECTOR_SIZE; i++)
    {
        if ((mask & (1 << i)) != 0x0)
        {
            c++;
        }
    }

    return c;
}

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
guessp_16(__m512 dl,
          __m512 ul,
          __m512 pl,
          __m512 cl,
          __m512 dr,
          __m512 ur,
          __m512 pr,
          __m512 cr,
          __m512 *pm)
{
    __m512 two, half, cup, ppv, pmin, pmax, qmax, pq, um, ptl, ptr, gel, ger, pqcr;
    __mmask16 cond_pvrs, cond_ppv, ncond_ppv;

    two = SET1(2.0);
    half = SET1(0.5);
    cup = MUL(SET1(0.25), MUL(ADD(dl, dr), ADD(cl, cr)));
    ppv = MUL(half, FMADD(SUB(ul, ur), cup, ADD(pl, pr)));
    ppv = MAX(ppv, z);
    pmin = MIN(pl, pr);
    pmax = MAX(pl, pr);
    qmax = DIV(pmax, pmin);

    // Conditions.
    cond_pvrs = CMP(qmax, two, _MM_CMPINT_LE)
                & CMP(pmin, ppv, _MM_CMPINT_LE)
                & CMP(ppv, pmax, _MM_CMPINT_LE);
    cond_ppv = _mm512_mask_cmp_ps_mask(~cond_pvrs, ppv, pmin, _MM_CMPINT_LT);
    ncond_ppv = ~cond_pvrs & ~cond_ppv;

    // The first branch.
    *pm = _mm512_mask_mov_ps(*pm, cond_pvrs, ppv);

    // The second branch.
    if (cond_ppv != 0x0)
    {
        pq = _mm512_mask_pow_ps(z, cond_ppv,
                                _mm512_mask_div_ps(z, cond_ppv, pl, pr), g1);
        pqcr = MUL(pq, cr);
        um = _mm512_mask_div_ps(z, cond_ppv,
                                FMADD(FMADD(SUB(pqcr, cr), g4, ur), cl, MUL(pqcr, ul)),
                                ADD(pqcr, cl));
        ptl = FMADD(_mm512_mask_div_ps(z, cond_ppv, SUB(ul, um), cl), g7, one);
        ptr = FMADD(_mm512_mask_div_ps(z, cond_ppv, SUB(um, ur), cr), g7, one);
        *pm = _mm512_mask_mul_ps(*pm, cond_ppv, half,
                                 ADD(_mm512_mask_pow_ps(z, cond_ppv, MUL(pl, ptl), g3),
                                     _mm512_mask_pow_ps(z, cond_ppv, MUL(pr, ptr), g3)));
    }

    // The third branch.
    if (ncond_ppv != 0x0)
    {
        gel = SQRT(_mm512_mask_div_ps(z, ncond_ppv, g5, MUL(FMADD(g6, pl, ppv), dl)));
        ger = SQRT(_mm512_mask_div_ps(z, ncond_ppv, g5, MUL(FMADD(g6, pr, ppv), dr)));
        *pm = _mm512_mask_div_ps(*pm, ncond_ppv,
                                 FMADD(gel, pl, FMADD(ger, pr, SUB(ul, ur))),
                                 ADD(gel, ger));
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
/// \param[in] m - mask for operations
static void
prefun_16(__m512 *f,
          __m512 *fd,
          __m512 p,
          __m512 dk,
          __m512 pk,
          __m512 ck,
          __mmask16 m)
{
    __m512 pratio, ak, bkp, ppk, qrt;
    __mmask16 cond, ncond;

    // Conditions.
    cond = _mm512_mask_cmp_ps_mask(m, p, pk, _MM_CMPINT_LE);
    ncond = m & ~cond;

    // The first branch.
    if (cond != 0x0)
    {
        pratio = _mm512_mask_div_ps(z, cond, p, pk);
        *f = _mm512_mask_mul_ps(*f, cond, MUL(g4, ck),
                                SUB(_mm512_mask_pow_ps(z, cond, pratio, g1), one));
        *fd = _mm512_mask_div_ps(*fd, cond,
                                 _mm512_mask_pow_ps(z, cond, pratio, SUB(z, g2)),
                                 MUL(dk, ck));
    }

    // The second branch.
    if (ncond != 0x0)
    {
        ak = _mm512_mask_div_ps(z, ncond, g5, dk);
        bkp = FMADD(g6, pk, p);
        ppk = SUB(p, pk);
        qrt = _mm512_mask_sqrt_ps(z, ncond,
                                  _mm512_mask_div_ps(z, ncond, ak, bkp));
        *f = _mm512_mask_mul_ps(*f, ncond, ppk, qrt);
        *fd = _mm512_mask_mul_ps(*fd, ncond, qrt,
                                 _mm512_fnmadd_ps(_mm512_mask_div_ps(z, ncond, ppk, bkp),
                                                  SET1(0.5), one));
    }
}

/// \brief Pressure and speed calculation in star region.
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
/// \param[out] p - pressure in star region
/// \param[out] u - velocity in star region
static void
starpu_16(__m512 dl,
          __m512 ul,
          __m512 pl,
          __m512 cl,
          __m512 dr,
          __m512 ur,
          __m512 pr,
          __m512 cr,
          __m512 *p,
          __m512 *u)
{
    __m512 two, tolpre, tolpre2, udiff, pold, fl, fld, fr, frd, change;
    __mmask16 cond_break, cond_neg, m;
    const int nriter = 20;
    int iter = 1;

    two = SET1(2.0);
    tolpre = SET1(1.0e-6);
    tolpre2 = SET1(5.0e-7);
    udiff = SUB(ur, ul);

    guessp_16(dl, ul, pl, cl, dr, ur, pr, cr, &pold);

    // Start with full mask.
    m = 0xFFFF;

    for (; (iter <= nriter) && (m != 0x0); iter++)
    {
        prefun_16(&fl, &fld, pold, dl, pl, cl, m);
        prefun_16(&fr, &frd, pold, dr, pr, cr, m);
        *p = _mm512_mask_sub_ps(*p, m, pold,
                                _mm512_mask_div_ps(z, m,
                                                   ADD(ADD(fl, fr), udiff),
                                                   ADD(fld, frd)));
        change = ABS(_mm512_mask_div_ps(z, m, SUB(*p, pold), ADD(*p, pold)));
        cond_break = _mm512_mask_cmp_ps_mask(m, change, tolpre2, _MM_CMPINT_LE);
        m &= ~cond_break;
        cond_neg = _mm512_mask_cmp_ps_mask(m, *p, z, _MM_CMPINT_LT);
        *p = _mm512_mask_mov_ps(*p, cond_neg, tolpre);
        pold = _mm512_mask_mov_ps(pold, m, *p);
    }

    // Check for divergence.
    if (iter > nriter)
    {
        cout << "divergence in Newton-Raphson iteration" << endl;
        exit(1);
    }

    *u = MUL(SET1(0.5), ADD(ADD(ul, ur), SUB(fr, fl)));
}

/// \brief Final analyze of the configuration.
///
/// Purpose is to sample the solution throughout the wave
/// pattern. Pressure pm and velocit
/// star region are known. Sampling is performed
/// in terms of the 'speed' s = x/t. Sampled
/// values are d, u, p.
///
/// \param[in] dl - left side density
/// \param[in] ul - left side velocity x component
/// \param[in] vl - left side velocity y component
/// \param[in] wl - left side velocity z component
/// \param[in] pl - left side pressure
/// \param[in] cl - left side sound velocity
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity x component
/// \param[in] vr - right side velocity y component
/// \param[in] wr - right side velocity z component
/// \param[in] pr - right side pressure
/// \param[in] cr - right side sound velocity
/// \param[out] d - result density
/// \param[out] u - result velocity x component
/// \param[out] v - result velocity y component
/// \param[out] w - result velocity z component
/// \param[out] p - result pressure
static void
sample_16(__m512 dl,
          __m512 ul,
          __m512 vl,
          __m512 wl,
          __m512 pl,
          __m512 cl,
          __m512 dr,
          __m512 ur,
          __m512 vr,
          __m512 wr,
          __m512 pr,
          __m512 cr,
          __m512 pm,
          __m512 um,
          __m512 *d,
          __m512 *u,
          __m512 *v,
          __m512 *w,
          __m512 *p)
{
    __m512 c, ums, pms, sh, st, s, uc;
    __mmask16 cond_um, cond_pm, cond_sh, cond_st, cond_s, cond_sh_st;

    // d/u/p/c/ums
    cond_um = _mm512_cmp_ps_mask(um, z, _MM_CMPINT_LT);
    *d = _mm512_mask_blend_ps(cond_um, dl, dr);
    *u = _mm512_mask_blend_ps(cond_um, ul, ur);
    *v = _mm512_mask_blend_ps(cond_um, vl, vr);
    *w = _mm512_mask_blend_ps(cond_um, wl, wr);
    *p = _mm512_mask_blend_ps(cond_um, pl, pr);
    c = _mm512_mask_blend_ps(cond_um, cl, cr);
    ums = um;
    *u = _mm512_mask_sub_ps(*u, cond_um, z, *u);
    ums = _mm512_mask_sub_ps(ums, cond_um, z, ums);

    // Calculate main values.
    pms = DIV(pm, *p);
    sh = SUB(*u, c);
    st = _mm512_fnmadd_ps(POW(pms, g1), c, ums);
    s = _mm512_fnmadd_ps(c, SQRT(FMADD(g2, pms, g1)), *u);

    // Conditions.
    cond_pm = _mm512_cmp_ps_mask(pm, *p, _MM_CMPINT_LE);
    cond_sh = _mm512_mask_cmp_ps_mask(cond_pm, sh, z, _MM_CMPINT_LT);
    cond_st = _mm512_mask_cmp_ps_mask(cond_sh, st, z, _MM_CMPINT_LT);
    cond_s = _mm512_mask_cmp_ps_mask(~cond_pm, s, z, _MM_CMPINT_LT);

    // Store.
    *d = _mm512_mask_mov_ps(*d, cond_st, MUL(*d, POW(pms, SET1(1.0 / GAMA))));
    *d = _mm512_mask_mov_ps(*d, cond_s, MUL(*d, DIV(ADD(pms, g6), FMADD(pms, g6, one))));
    *u = _mm512_mask_mov_ps(*u, cond_st | cond_s, ums);
    *p = _mm512_mask_mov_ps(*p, cond_st | cond_s, pm);

    // Low prob - ignnore it.
    cond_sh_st = cond_sh & ~cond_st;
    if (cond_sh_st != 0x0)
    {
        *u = _mm512_mask_mov_ps(*u, cond_sh_st, MUL(g5, FMADD(g7, *u, c)));
        uc = DIV(*u, c);
        *d = _mm512_mask_mov_ps(*d, cond_sh_st, MUL(*d, POW(uc, g4)));
        *p = _mm512_mask_mov_ps(*p, cond_sh_st, MUL(*p, POW(uc, g3)));
    }

    // Final store.
    *u = _mm512_mask_sub_ps(*u, cond_um, z, *u);
}

/// \brief Riemann solver for 16 cases.
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
/// \param[out] v - result velocity reference x component
/// \param[out] w - result velocity reference x component
/// \param[out] p - result pressure reference
static void
riemann_16(__m512 dl,
           __m512 ul,
           __m512 vl,
           __m512 wl,
           __m512 pl,
           __m512 dr,
           __m512 ur,
           __m512 vr,
           __m512 wr,
           __m512 pr,
           __m512 *d,
           __m512 *u,
           __m512 *v,
           __m512 *w,
           __m512 *p)
{
    __m512 cl, cr, pm, um;
    __mmask16 vacuum_mask;

    cl = SQRT(DIV(MUL(gama, pl), dl));
    cr = SQRT(DIV(MUL(gama, pr), dr));
    vacuum_mask = CMP(MUL(g4, ADD(cl, cr)),
                      SUB(ur, ul),
                      _MM_CMPINT_LE);

    // Vacuum check.
    if (vacuum_mask != 0x0)
    {
        cerr << "VACUUM" << endl;
        exit(1);
    }

    starpu_16(dl, ul, pl, cl, dr, ur, pr, cr, &pm, &um);
    sample_16(dl, ul, vl, wl, pl, cl,
              dr, ur, vr, wr, pr, cr,
              pm, um,
              d, u, v, w, p);
}

#endif // INTEL

/// \brief Riemann solver.
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
void
riemann_opt(int c,
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

#ifndef INTEL

    riemann(c,
            dl, ul, vl, wl, pl,
            dr, ur, vr, wr, pr,
            d, u, v, w, p,
            nt);

#else

    __assume_aligned(dl, 64);
    __assume_aligned(ul, 64);
    __assume_aligned(vl, 64);
    __assume_aligned(wl, 64);
    __assume_aligned(pl, 64);
    __assume_aligned(dr, 64);
    __assume_aligned(ur, 64);
    __assume_aligned(vr, 64);
    __assume_aligned(wr, 64);
    __assume_aligned(pr, 64);
    __assume_aligned(d, 64);
    __assume_aligned(u, 64);
    __assume_aligned(v, 64);
    __assume_aligned(w, 64);
    __assume_aligned(p, 64);

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
        __m512 vd, vu, vv, vw, vp;

        for (int i = lb * FP16_VECTOR_SIZE;
             i < ub * FP16_VECTOR_SIZE;
             i += FP16_VECTOR_SIZE)
        {
            riemann_16(LD(dl + i), LD(ul + i), LD(vl + i), LD(wl + i), LD(pl + i),
                       LD(dr + i), LD(ur + i), LD(vr + i), LD(wer + i), LD(pr + i),
                       &vd, &vu, &vv, &vw, &vp);
            ST(d + i, vd);
            ST(u + i, vu);
            ST(v + i, vv);
            ST(w + i, vw);
            ST(p + i, vp);
        }
    }

    omp_set_num_threads(1);

    // Tail.
    riemann(c_tail,
            dl + c_base, ul + c_base, vl + c_base, wl + c_base, pl + c_base,
            dr + c_base, ur + c_base, vr + c_base, wr + c_base, pr + c_base,
            d + c_base, u + c_base, v + c_base, w + c_base, p + c_base);

#endif // INTEL

}
