/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_H
#define RIEMANN_H

void init_gamas();
void riemann(int c,
             float *dl, float *ul, float *pl,
             float *dr, float *ur, float *pr,
             float *d, float *u, float *p);

#endif
