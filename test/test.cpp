/// \file
/// \brief Test.
///
/// Test for riemann.

#include "riemann.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <omp.h>
using namespace std;

/// \brief Test data <c>dl</c>.
double dls[] =
{

#include "test_data_dl.inc"

};

/// \brief Test data <c>ul</c>.
double uls[] =
{

#include "test_data_ul.inc"

};

/// \brief Test data <c>pl</c>.
double pls[] =
{

#include "test_data_pl.inc"

};

/// \brief Test data <c>dr</c>.
double drs[] =
{

#include "test_data_dr.inc"

};

/// \brief Test data <c>ur</c>.
double urs[] =
{

#include "test_data_ur.inc"

};

/// \brief Test data <c>pr</c>.
double prs[] =
{

#include "test_data_pr.inc"

};

/// \brief Test data <c>d</c>.
double ds[] =
{

#include "test_data_d.inc"

};

/// \brief Test data <c>u</c>.
double us[] =
{

#include "test_data_u.inc"

};

/// \brief Test data <c>p</c>.
double ps[] =
{

#include "test_data_p.inc"

};

/// \brief Test.
int main()
{
    int test_cases = sizeof(dls) / (sizeof(dls[0]));
    double d, u, p;
    double e = 1e-4;

    init_gamas();

    cout << "test begin : " << test_cases << " test cases" << endl;

    double t_start = omp_get_wtime();

    for (int i = 0; i < test_cases; i++)
    {
        riemann(dls[i], uls[i], pls[i],
                drs[i], urs[i], prs[i],
                d, u, p);

        double diff_d = abs(d - ds[i]);
        double diff_u = abs(u - us[i]);
        double diff_p = abs(p - ps[i]);

        if (!(diff_d < e) && (diff_u < e) && (diff_p < e))
        {
            cerr << "error : " << endl;
            cerr << "  res : " << d << ", " << u << ", " << p << endl;
            cerr << "right : " << ds[i] << ", " << us[i] << ", " << ps[i] << endl; 
            cerr << " diff : " << diff_d << ", " << diff_u << ", " << diff_p << endl;
            exit(1);
        }
    }

    double t_end = omp_get_wtime();
    double t_len = t_end - t_start;

    cout << "test done : " << (t_end - t_start) << " seconds" << endl;

    return 0;
}
