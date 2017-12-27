/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_OPT_H
#define RIEMANN_OPT_H

void init_gamas_opt();
void riemann_opt(int c,
                 float *dl, float *ul, float *pl,
                 float *dr, float *ur, float *pr,
                 float *d, float *u, float *p);

#endif
