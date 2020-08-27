/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_H
#define RIEMANN_H

/// \brief Gama value.
#define GAMA 1.4

/// \brief Gama special constant 1.
#define G1 (((GAMA) - 1.0) / (2.0 * (GAMA)))

/// \brief Gama special constant 2.
#define G2 (((GAMA) + 1.0) / (2.0 * (GAMA)))

/// \brief Gama special constant 3.
#define G3 (2.0 * (GAMA) / ((GAMA) - 1.0))

/// \brief Gama special constant 4.
#define G4 (2.0 / ((GAMA) - 1.0))

/// \brief Gama special constant 5.
#define G5 (2.0 / ((GAMA) + 1.0))

/// \brief Gama special constant 6.
#define G6 (((GAMA) - 1.0) / ((GAMA) + 1.0))

/// \brief Gama special constant 7.
#define G7 (((GAMA) - 1.0) / 2.0)

/// \brief Gama special constant 8.
#define G8 ((GAMA) - 1.0)

//
// Prototypes.
//

// Single not optimized original version.
void
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
        float &p);

// Not optimized version for multiple data.
void
riemann(int c,
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
        float *p);

// Vectorized version for miltiple data.
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
            float *p);

#endif // !RIEMANN_H
