/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_H
#define RIEMANN_H

/// \brief Gama value.
#define GAMA 1.4

/// \brief Gama 1.
#define G1 (((GAMA) - 1.0) / (2.0 * (GAMA)))

/// \brief Gama 2.
#define G2 (((GAMA) + 1.0) / (2.0 * (GAMA)))

/// \brief Gama 3.
#define G3 (2.0 * (GAMA) / ((GAMA) - 1.0))

/// \brief Gama 4.
#define G4 (2.0 / ((GAMA) - 1.0))

/// \brief Gama 5.
#define G5 (2.0 / ((GAMA) + 1.0))

/// \brief Gama 6.
#define G6 (((GAMA) - 1.0) / ((GAMA) + 1.0))

/// \brief Gama 7.
#define G7 (((GAMA) - 1.0) / 2.0)

/// \brief Gama 8.
#define G8 ((GAMA) - 1.0)

// Prototypes.
void riemann(int c,
             float *dl, float *ul, float *pl,
             float *dr, float *ur, float *pr,
             float *d, float *u, float *p);
void riemann_opt(int c,
                 float *dl, float *ul, float *pl,
                 float *dr, float *ur, float *pr,
                 float *d, float *u, float *p);

#endif
