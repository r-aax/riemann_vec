/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_OPT_H
#define RIEMANN_OPT_H

void init_gamas_opt();
void riemann_opt(int c,
                 double *dl, double *ul, double *pl,
                 double *dr, double *ur, double *pr,
                 double *d, double *u, double *p);

#endif
