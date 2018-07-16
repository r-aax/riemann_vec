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
#include <iomanip>
#include <omp.h>
using namespace std;

/// \brief Repeats count.
#define REPEATS 5

/// \brief Test data <c>dl</c>.
float dls[] =
{

#include "test_data_dl.inc"

};

/// \brief Test data <c>ul</c>.
float uls[] =
{

#include "test_data_ul.inc"

};

/// \brief Test data <c>pl</c>.
float pls[] =
{

#include "test_data_pl.inc"

};

/// \brief Test data <c>dr</c>.
float drs[] =
{

#include "test_data_dr.inc"

};

/// \brief Test data <c>ur</c>.
float urs[] =
{

#include "test_data_ur.inc"

};

/// \brief Test data <c>pr</c>.
float prs[] =
{

#include "test_data_pr.inc"

};

/// \brief Test data <c>d</c>.
float ds_orig[] =
{

#include "test_data_d.inc"

};

/// \brief Test data <c>u</c>.
float us_orig[] =
{

#include "test_data_u.inc"

};

/// \brief Test data <c>p</c>.
float ps_orig[] =
{

#include "test_data_p.inc"

};

/// \brief Calculated result <c>d</c>.
float *ds;

/// \brief Calculated result <c>u</c>.
float *us;

/// \brief Calculated result <c>p</c>.
float *ps;

/// \brief Test cases count.
int test_cases = sizeof(dls) / sizeof(dls[0]);

/// \brief Check function.
void check()
{
    float e = 1e-3;

    for (int i = 0; i < test_cases; i++)
    {
        float diff_d = abs(ds[i] - ds_orig[i]);
        float diff_u = abs(us[i] - us_orig[i]);
        float diff_p = abs(ps[i] - ps_orig[i]);

        if (!((diff_d < e) && (diff_u < e) && (diff_p < e)))
        {
            cerr << "error : " << endl;
            cerr << "  res : " << ds[i] << ", " << us[i] << ", " << ps[i] << endl;
            cerr << "right : " << ds_orig[i] << ", " << us_orig[i] << ", " << ps_orig[i] << endl; 
            cerr << " diff : " << diff_d << ", " << diff_u << ", " << diff_p << endl;
            exit(1);
        }
    }
}

/// \brief Clean result data.
void clean()
{
    for (int i = 0; i < test_cases; i++)
    {
        ds[i] = 0.0;
        us[i] = 0.0;
        ps[i] = 0.0;
    }
}

/// \brief Run riemann solver test and print information.
///
/// \param[in] solver - solver function
/// \param[in] str - description
///
/// \return
/// Execution time.
double run(void (*solver)(int,
                          float *, float *, float *,
                          float *, float *, float *,
                          float *, float *, float *),
           string str)
{

/// \brief Inner repeats count.
#define INNER_REPEATS 10

    clean();
    double t_start = omp_get_wtime();
    for (int i = 0; i < INNER_REPEATS; i++)
    {
        solver(test_cases, dls, uls, pls, drs, urs, prs, ds, us, ps);
    }
    double t_end = omp_get_wtime();
    check();    
    double t_len = t_end - t_start;
    cout << setw(15) << str << " ~ " << t_len << " seconds" << endl;

    return t_len;
}

/// \brief Min value in array.
///
/// \param d - array
/// \param c - element count
///
/// \return
/// Min value.
double array_min(double *d, int c)
{
    double m = d[0];

    for (int i = 1; i < c; i++)
    {
        if (d[i] < m)
        {
            m = d[i];
        }
    }

    return m;
}

/// \brief Test.
int main()
{
    double times[REPEATS];
    double times_opt[REPEATS];

    if (!((sizeof(dls) == sizeof(uls)) && (sizeof(uls) == sizeof(pls))
          && (sizeof(pls) == sizeof(drs)) && (sizeof(drs) == sizeof(urs))
          && (sizeof(urs) == sizeof(prs)) && (sizeof(prs) == sizeof(ds_orig))
          && (sizeof(ds_orig) == sizeof(us_orig)) && (sizeof(us_orig) == sizeof(ps_orig))))
    {
        cout << "error : test data corrupted" << endl;
        exit(1);
    }

    ds = new float[test_cases];
    us = new float[test_cases];
    ps = new float[test_cases];

    cout << "test begin : " << test_cases << " test cases" << endl;

    for (int i = 0; i < REPEATS; i++)
    {
        times[i] = run(riemann, "not optimized");
    }

    cout << "----------" << endl;

    for (int i = 0; i < REPEATS; i++)
    {
        times_opt[i] = run(riemann_opt, "optimized");
    }

    double min_time = array_min(times, REPEATS);
    double min_time_opt = array_min(times_opt, REPEATS);
    double time_reduce = ((min_time - min_time_opt) / min_time) * 100.0;
    double speedup_x = min_time / min_time_opt;
    cout << "test done : time_reduce = " << setprecision(2) << time_reduce
         << "%, speedup_x = " << setprecision(3) << speedup_x << endl;

    // Free memory.
    delete ds;
    delete us;
    delete ps;

    return 0;
}
