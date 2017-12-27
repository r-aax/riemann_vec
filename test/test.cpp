/// \file
/// \brief Test.
///
/// Test for riemann.

#include "riemann.h"
#include "riemann_opt.h"
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
    float e = 1e-4;

    for (int i = 0; i < test_cases; i++)
    {
        float diff_d = abs(ds[i] - ds_orig[i]);
        float diff_u = abs(us[i] - us_orig[i]);
        float diff_p = abs(ps[i] - ps_orig[i]);

        if (!(diff_d < e) && (diff_u < e) && (diff_p < e))
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
/// \param[in] init - init function
/// \param[in] solver - solver function
/// \param[in] str - description
void run(void (*init)(),
         void (*solver)(int,
                        float *, float *, float *,
                        float *, float *, float *,
                        float *, float *, float *),
         string str)
{
    clean();
    init();
    double t_start = omp_get_wtime();
    solver(test_cases, dls, uls, pls, drs, urs, prs, ds, us, ps);
    double t_end = omp_get_wtime();
    check();    
    double t_len = t_end - t_start;
    cout << setw(15) << str << " ~ " << (t_end - t_start) << " seconds" << endl;
}

/// \brief Test.
int main()
{
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

    init_gamas();

    cout << "test begin : " << test_cases << " test cases" << endl;

    for (int i = 0; i < REPEATS; i++)
    {
        run(init_gamas, riemann, "not optimized");
    }

    cout << "----------" << endl;

    for (int i = 0; i < REPEATS; i++)
    {
        run(init_gamas_opt, riemann_opt, "optimized");
    }

    cout << "test done" << endl;

    // Free memory.
    delete ds;
    delete us;
    delete ps;

    return 0;
}
