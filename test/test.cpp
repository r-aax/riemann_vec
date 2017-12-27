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
double ds_orig[] =
{

#include "test_data_d.inc"

};

/// \brief Test data <c>u</c>.
double us_orig[] =
{

#include "test_data_u.inc"

};

/// \brief Test data <c>p</c>.
double ps_orig[] =
{

#include "test_data_p.inc"

};

/// \brief Test.
int main()
{
    int dls_count = sizeof(dls) / (sizeof(dls[0]));
    int uls_count = sizeof(uls) / (sizeof(uls[0]));
    int pls_count = sizeof(pls) / (sizeof(pls[0]));
    int drs_count = sizeof(drs) / (sizeof(drs[0]));
    int urs_count = sizeof(urs) / (sizeof(urs[0]));
    int prs_count = sizeof(prs) / (sizeof(prs[0]));
    int ds_count = sizeof(ds_orig) / (sizeof(ds_orig[0]));
    int us_count = sizeof(us_orig) / (sizeof(us_orig[0]));
    int ps_count = sizeof(ps_orig) / (sizeof(ps_orig[0]));

    if (!((dls_count == uls_count) && (uls_count == pls_count)
          && (pls_count == drs_count) && (drs_count == urs_count)
          && (urs_count == prs_count) && (prs_count == ds_count)
          && (ds_count == us_count) && (us_count == ps_count)))
    {
        cout << "error : test data corrupted" << endl;
        exit(1);
    }

    int test_cases = dls_count;
    double *ds, *us, *ps;
    double e = 1e-4;

    ds = new double[test_cases];
    us = new double[test_cases];
    ps = new double[test_cases];

    init_gamas();

    cout << "test begin : " << test_cases << " test cases" << endl;

    double t_start = omp_get_wtime();

    // Calculations loop.
    for (int i = 0; i < test_cases; i++)
    {
        double d, u, p;

        riemann(dls[i], uls[i], pls[i],
                drs[i], urs[i], prs[i],
                d, u, p);

        ds[i] = d;
        us[i] = u;
        ps[i] = p;
    }

    double t_end = omp_get_wtime();

    // Check loop.
    for (int i = 0; i < test_cases; i++)
    {
        double diff_d = abs(ds[i] - ds_orig[i]);
        double diff_u = abs(us[i] - us_orig[i]);
        double diff_p = abs(ps[i] - ps_orig[i]);

        if (!(diff_d < e) && (diff_u < e) && (diff_p < e))
        {
            cerr << "error : " << endl;
            cerr << "  res : " << ds[i] << ", " << us[i] << ", " << ps[i] << endl;
            cerr << "right : " << ds_orig[i] << ", " << us_orig[i] << ", " << ps_orig[i] << endl; 
            cerr << " diff : " << diff_d << ", " << diff_u << ", " << diff_p << endl;
            exit(1);
        }
    }

    double t_len = t_end - t_start;
    cout << "test done : " << (t_end - t_start) << " seconds" << endl;

    return 0;
}
