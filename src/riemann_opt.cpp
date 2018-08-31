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

/// \brief Short name for load intrinsic.
#define LD(ADDR) _mm512_load_ps(ADDR)

/// \brief Short name for store intrinsic.
#define ST(ADDR, VAL) _mm512_store_ps(ADDR, VAL)

/// \brief Short name for add intrinsic.
#define ADD(va, vb) _mm512_add_ps(va, vb)

/// \brief Short name for mul intrinsic.
#define MUL(va, vb) _mm512_mul_ps(va, vb)

/// \brief Short name for sub intrinsic.
#define SUB(va, vb) _mm512_sub_ps(va, vb)

/// \brief Short name for div intrinsic.
#define DIV(va, vb) _mm512_div_ps(va, vb)

/// \brief Short name for pow intrinsic.
#define POW(va, vb) _mm512_pow_ps(va, vb)

/// \brief Short name for sqrt intrinsic.
#define SQRT(va) _mm512_sqrt_ps(va)

/// \brief Short name for abs intrinsic.
#define ABS(v) _mm512_abs_ps(v)

/// \brief Short name for max intrinsic.
#define MAX(va, vb) _mm512_max_ps(va, vb)

/// \brief Short name for min intrinsic.
#define MIN(va, vb) _mm512_min_ps(va, vb)

/// \bries Short name for cmp.
#define CMP(va, vb, cmp) _mm512_cmp_ps_mask(va, vb, cmp)

/// \brief Short name for setzero intrinsic.
#define SETZERO() _mm512_setzero_ps()

/// \brief Short name for set1 intrinsic.
#define SET1(v) _mm512_set1_ps(v)

/// \brief Zero.
__m512 z = SETZERO();

/// \brief 1.
__m512 v1 = SET1(1.0);

/// \brief GAMA.
__m512 gama = SET1(GAMA);

/// \brief 1/GAMA.
__m512 igama = SET1(1.0 / GAMA);

/// \brief G1.
__m512 g1 = SET1(G1);

/// \brief G2.
__m512 g2 = SET1(G2);

/// \brief G3.
__m512 g3 = SET1(G3);

/// \brief G4.
__m512 g4 = SET1(G4);

/// \brief G5.
__m512 g5 = SET1(G5);

/// \brief G6.
__m512 g6 = SET1(G6);

/// \brief G7.
__m512 g7 = SET1(G7);

/// \brief Constant for starpu.
__m512 tolpre = SET1(1.0e-6);

/// \brief Get <c>i</c>-th element from vector.
///
/// \param[in] v - vector
/// \param[in] i - index
///
/// \return
/// Element.
static float Get(__m512 v, int i)
{
    float arr[16];

    ST(&arr[0], v);

    return arr[i];
}

/// \brief Set <c>i</c>-th element in vector.
///
/// \param[in,out] v - vector
/// \param[in] i - index
/// \param[in] f - value
static void Set(__m512 *v, int i, float f)
{
    float arr[16];

    ST(&arr[0], *v);
    arr[i] = f;
    *v = LD(&arr[0]);
}

/// \brief Print.
///
/// \param[in] v - vector
static void Print(__m512 v)
{
    float arr[16];

    ST(&arr[0], v);

    printf("[");
    for (int i = 0; i < 16; i++)
    {
        printf(" %f ", arr[i]);
    }
    printf("]\n");
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
static void guessp_16(__m512 dl, __m512 ul, __m512 pl, __m512 cl,
                      __m512 dr, __m512 ur, __m512 pr, __m512 cr,
                      __m512 *pm)
{
    __m512 half = SET1(0.5);
    __m512 quser = SET1(2.0);
    __m512 cup, ppv, pmin, pmax, qmax, pq, um, ptl, ptr, gel, ger;
    __mmask16 cond_pvrs, cond_ppv, ncond_ppv;

    cup = MUL(half, MUL(ADD(dl, dr), ADD(cl, cr)));
    ppv = MUL(half, ADD(ADD(pl, pr), MUL(cup, SUB(ul, ur))));
    ppv = MAX(ppv, z);
    pmin = MIN(pl, pr);
    pmax = MAX(pl, pr);
    qmax = DIV(pmax, pmin);
    cond_pvrs = CMP(qmax, quser, _MM_CMPINT_LE)
                && CMP(pmin, ppv, _MM_CMPINT_LE)
                && CMP(ppv, pmax, _MM_CMPINT_LE);
    *pm = _mm512_mask_mov_ps(*pm, cond_pvrs, ppv);
    cond_ppv = _mm512_mask_cmp_ps_mask(~cond_pvrs, ppv, pmin, _MM_CMPINT_LT);
    ncond_ppv = ~cond_pvrs & ~cond_ppv;
    pq = _mm512_mask_pow_ps(z, cond_ppv,
                            _mm512_mask_div_ps(z, cond_ppv, pl, pr), g1);
    um = _mm512_mask_div_ps(z, cond_ppv,
                            ADD(ADD(MUL(pq, _mm512_mask_div_ps(z, cond_ppv, ul, cl)),
                                    _mm512_mask_div_ps(z, cond_ppv, ur, cr)),
                                MUL(g4, SUB(pq, v1))),
                            ADD(_mm512_mask_div_ps(z, cond_ppv, pq, cl),
                                _mm512_mask_div_ps(z, cond_ppv, v1, cr)));
    ptl = ADD(v1, MUL(g7, _mm512_mask_div_ps(z, cond_ppv, SUB(ul, um), cl)));
    ptr = ADD(v1, MUL(g7, _mm512_mask_div_ps(z, cond_ppv, SUB(um, ur), cr)));
    *pm = _mm512_mask_mul_ps(*pm, cond_ppv, half,
                             ADD(_mm512_mask_pow_ps(z, cond_ppv, MUL(pl, ptl), g3),
                                 _mm512_mask_pow_ps(z, cond_ppv, MUL(pr, ptr), g3)));
    gel = _mm512_mask_sqrt_ps(z, ncond_ppv,
                              _mm512_mask_div_ps(z, ncond_ppv,
                                                 _mm512_mask_div_ps(z, ncond_ppv, g5, dl),
                                                 ADD(MUL(g6, pl), ppv)));
    ger = _mm512_mask_sqrt_ps(z, ncond_ppv,
                              _mm512_mask_div_ps(z, ncond_ppv,
                                                 _mm512_mask_div_ps(z, ncond_ppv, g5, dr),
                                                 ADD(MUL(g6, pr), ppv)));
    *pm = _mm512_mask_div_ps(*pm, ncond_ppv,
                             SUB(ADD(MUL(gel, pl), MUL(ger, pr)), SUB(ur, ul)),
                             ADD(gel, ger));
}

/// \brief
///
/// Purpose is to evaluate the pressure functions
/// fl and fr in exact Riemann solver
/// and their first derivatives.
///
/// \param[in,out] f - ?
/// \param[in,out] fd - ?
/// \param[in] p - ?
/// \param[in] dk - ?
/// \param[in] pk - ?
/// \param[in] ck - ?
/// \param[in] m - mask for operations
static void prefun_16(__m512 *f, __m512 *fd, __m512 p,
                      __m512 dk, __m512 pk, __m512 ck,
                      __mmask16 m)
{
    __mmask16 cond = _mm512_mask_cmp_ps_mask(m, p, pk, _MM_CMPINT_LE);
    __mmask16 ncond = m & ~cond;
    __m512 pratio = _mm512_mask_div_ps(z, cond, p, pk);
    __m512 ak = _mm512_mask_div_ps(z, ncond, g5, dk);
    __m512 bk = MUL(g6, pk);
    __m512 qrt = _mm512_mask_sqrt_ps(z, ncond,
                                     _mm512_mask_div_ps(z, ncond, ak, ADD(bk, p)));

    // Output values.
    *f = _mm512_mask_mul_ps(*f, cond, MUL(g4, ck),
                            SUB(_mm512_mask_pow_ps(z, cond, pratio, g1), v1));
    *fd = _mm512_mask_mul_ps(*fd, cond,
                             _mm512_mask_div_ps(z, cond, v1, MUL(dk, ck)),
                             _mm512_mask_pow_ps(z, cond, pratio, SUB(z, g2)));
    *f = _mm512_mask_mul_ps(*f, ncond, SUB(p, pk), qrt);
    *fd = _mm512_mask_mul_ps(*fd, ncond, qrt,
                             SUB(v1, MUL(SET1(0.5),
                                         _mm512_mask_div_ps(z, ncond,
                                                            SUB(p, pk),
                                                            ADD(bk, p)))));
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
/// \param[out] p - pressure in star region
/// \param[out] u - velocity in star region
static void starpu_16(__m512 dl, __m512 ul, __m512 pl, __m512 cl,
                      __m512 dr, __m512 ur, __m512 pr, __m512 cr,
                      __m512 *p, __m512 *u)
{
    __m512 pold, change, fl, fld, fr, frd;
    __m512 udiff = SUB(ur, ul);
    __mmask16 cond_break, cond_neg;
    __mmask16 m = 0xFFFF;
    const int nriter = 20;
    int iter = 1;

    guessp_16(dl, ul, pl, cl, dr, ur, pr, cr, &pold);

    for (; (iter <= nriter) && (m != 0x0); iter++)
    {
        prefun_16(&fl, &fld, pold, dl, pl, cl, m);
        prefun_16(&fr, &frd, pold, dr, pr, cr, m);
        *p = _mm512_mask_sub_ps(*p, m, pold,
                                _mm512_mask_div_ps(z, m,
                                                   ADD(ADD(fl, fr), udiff),
                                                   ADD(fld, frd)));
        change = MUL(SET1(2.0),
                     ABS(_mm512_mask_div_ps(z, m,
                                            SUB(*p, pold),
                                            ADD(*p, pold))));
        cond_break = _mm512_mask_cmp_ps_mask(m, change, tolpre, _MM_CMPINT_LE);
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
/// \param[out] d - result density
/// \param[out] u - result velocity
/// \param[out] p - result pressure
static void sample_16(__m512 dl, __m512 ul, __m512 pl, __m512 cl,
                      __m512 dr, __m512 ur, __m512 pr, __m512 cr,
                      __m512 pm, __m512 um,
                      __m512 *d, __m512 *u, __m512 *p)
{
    __m512 c, ums, pms, sh, st, s, uc;
    __mmask16 cond_um, cond_pm, cond_sh, cond_st, cond_s, cond_sh_st;

    // d/u/p/c/ums
    cond_um = _mm512_cmp_ps_mask(um, z, _MM_CMPINT_LT);
    *d = _mm512_mask_blend_ps(cond_um, dl, dr);
    *u = _mm512_mask_blend_ps(cond_um, ul, ur);
    *p = _mm512_mask_blend_ps(cond_um, pl, pr);
    c = _mm512_mask_blend_ps(cond_um, cl, cr);
    ums = um;
    *u = _mm512_mask_sub_ps(*u, cond_um, z, *u);
    ums = _mm512_mask_sub_ps(ums, cond_um, z, ums);

    // Calculate main values.
    pms = DIV(pm, *p);
    sh = SUB(*u, c);
    st = _mm512_fnmadd_ps(POW(pms, g1), c, ums);
    s = _mm512_fnmadd_ps(c, SQRT(_mm512_fmadd_ps(g2, pms, g1)), *u);

    // Conditions.
    cond_pm = _mm512_cmp_ps_mask(pm, *p, _MM_CMPINT_LE);
    cond_sh = _mm512_mask_cmp_ps_mask(cond_pm, sh, z, _MM_CMPINT_LT);
    cond_st = _mm512_mask_cmp_ps_mask(cond_sh, st, z, _MM_CMPINT_LT);
    cond_s = _mm512_mask_cmp_ps_mask(~cond_pm, s, z, _MM_CMPINT_LT);

    // Store.
    *d = _mm512_mask_mov_ps(*d, cond_st, MUL(*d, POW(pms, igama)));
    *d = _mm512_mask_mov_ps(*d, cond_s, MUL(*d, DIV(ADD(pms, g6), _mm512_fmadd_ps(pms, g6, v1))));
    *u = _mm512_mask_mov_ps(*u, cond_st | cond_s, ums);
    *p = _mm512_mask_mov_ps(*p, cond_st | cond_s, pm);

    // Low prob - ignnore it.
    cond_sh_st = cond_sh & ~cond_st;
    if (cond_sh_st)
    {
        *u = _mm512_mask_mov_ps(*u, cond_sh_st, MUL(g5, _mm512_fmadd_ps(g7, *u, c)));
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
/// \param[in] ul - left side velocity
/// \param[in] pl - left  side pressure
/// \param[in] dr - right side density
/// \param[in] ur - right side velocity
/// \param[in] pr - right side pressure
/// \param[out] d - result density reference
/// \param[out] u - result velocity reference
/// \param[out] p - result pressure reference
static void riemann_16(__m512 dl, __m512 ul, __m512 pl,
                       __m512 dr, __m512 ur, __m512 pr,
                       __m512 *d, __m512 *u, __m512 *p)
{
    __m512 cl = SQRT(DIV(MUL(gama, pl), dl));
    __m512 cr = SQRT(DIV(MUL(gama, pr), dr));
    __m512 pm, um;
    __mmask16 vacuum_mask = CMP(MUL(g4, ADD(cl, cr)),
                                SUB(ur, ul),
                                _MM_CMPINT_LE);

    // Vacuum check.
    if (vacuum_mask != 0x0)
    {
        cerr << "VACUUM" << endl;
        exit(1);
    }

    starpu_16(dl, ul, pl, cl, dr, ur, pr, cr, &pm, &um);
    sample_16(dl, ul, pl, cl, dr, ur, pr, cr, pm, um, d, u, p);
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

    __assume_aligned(dl, 64);
    __assume_aligned(ul, 64);
    __assume_aligned(pl, 64);
    __assume_aligned(dr, 64);
    __assume_aligned(ur, 64);
    __assume_aligned(pr, 64);
    __assume_aligned(d, 64);
    __assume_aligned(u, 64);
    __assume_aligned(p, 64);

    __m512 vd, vu, vp;
    int c_tail = c & 0xF;
    int c_base = c - c_tail;

    // Main body.
    for (int i = 0; i < c_base; i += 16)
    {
        riemann_16(LD(dl + i), LD(ul + i), LD(pl + i),
                   LD(dr + i), LD(ur + i), LD(pr + i),
                   &vd, &vu, &vp);
        ST(d + i, vd);
        ST(u + i, vu);
        ST(p + i, vp);
    }

    // Tail.
    riemann(c_tail,
            dl + c_base, ul + c_base, pl + c_base,
            dr + c_base, ur + c_base, pr + c_base,
            d + c_base, u + c_base, p + c_base);

#endif

}
