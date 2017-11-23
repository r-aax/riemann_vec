/// \file
/// \brief Functions declarations for Riemann solver.

#ifndef RIEMANN_H
#define RIEMANN_H

void init_gamas();
void riemann(double dl, double ul, double pl,
             double dr, double ur, double pr,
             double &d, double &u, double &p);

#endif
