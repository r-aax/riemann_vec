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

/// \brief Alignment of 64 bytes.
#ifdef INTEL
#define ALIGN_64 __declspec(align(64))
#else
#define ALIGN_64
#endif

/// \brief Test cases in big mode.
#define TEST_CASES_BIG 419996

/// \brief Test cases in small mode.
#define TEST_CASES_SMALL 32

#ifdef TEST_MODE
#if TEST_MODE == 0
#define TEST_CASES TEST_CASES_SMALL
#elif TEST_MODE == 1
#define TEST_CASES TEST_CASES_BIG
#else
#error "wrong TEST_MODE"
#endif
#else
#error "unknown TEST_MODE"
#endif

/// \brief Array length.
///
/// \param A - array
#define ARR_LEN(A) (sizeof(A) / sizeof(A[0]))

/// \brief Test data <c>dl</c>.
ALIGN_64 float dls[] =
{

#if TEST_MODE == 0
#include "small/test_data_dl.inc"
#else
#include "big/test_data_dl.inc"
#endif

};

/// \brief Test data <c>ul</c>.
ALIGN_64 float uls[] =
{

#if TEST_MODE == 0
#include "small/test_data_ul.inc"
#else
#include "big/test_data_ul.inc"
#endif

};

/// \brief Test data <c>vl</c>.
ALIGN_64 float vls[TEST_CASES] = { 0.0 };

/// \brief Test data <c>vl</c>.
ALIGN_64 float wls[TEST_CASES] = { 0.0 };

/// \brief Test data <c>pl</c>.
ALIGN_64 float pls[] =
{

#if TEST_MODE == 0
#include "small/test_data_pl.inc"
#else
#include "big/test_data_pl.inc"
#endif

};

/// \brief Test data <c>dr</c>.
ALIGN_64 float drs[] =
{

#if TEST_MODE == 0
#include "small/test_data_dr.inc"
#else
#include "big/test_data_dr.inc"
#endif

};

/// \brief Test data <c>ur</c>.
ALIGN_64 float urs[] =
{

#if TEST_MODE == 0
#include "small/test_data_ur.inc"
#else
#include "big/test_data_ur.inc"
#endif

};

/// \brief Test data <c>vr</c>.
ALIGN_64 float vrs[TEST_CASES] = { 0.0 };

/// \brief Test data <c>wr</c>.
ALIGN_64 float wrs[TEST_CASES] = { 0.0 };

/// \brief Test data <c>pr</c>.
ALIGN_64 float prs[] =
{

#if TEST_MODE == 0
#include "small/test_data_pr.inc"
#else
#include "big/test_data_pr.inc"
#endif

};

/// \brief Test data <c>d</c>.
ALIGN_64 float ds_orig[] =
{

#if TEST_MODE == 0
#include "small/test_data_d.inc"
#else
#include "big/test_data_d.inc"
#endif

};

/// \brief Test data <c>u</c>.
ALIGN_64 float us_orig[] =
{

#if TEST_MODE == 0
#include "small/test_data_u.inc"
#else
#include "big/test_data_u.inc"
#endif

};

/// \brief Test data <c>p</c>.
ALIGN_64 float ps_orig[] =
{

#if TEST_MODE == 0
#include "small/test_data_p.inc"
#else
#include "big/test_data_p.inc"
#endif

};

/// \brief Calculated result <c>d</c>.
ALIGN_64 float ds[TEST_CASES];

/// \brief Calculated result <c>u</c>.
ALIGN_64 float us[TEST_CASES];

/// \brief Calculated result <c>v</c>.
ALIGN_64 float vs[TEST_CASES];

/// \brief Calculated result <c>w</c>.
ALIGN_64 float ws_[TEST_CASES];

/// \brief Calculated result <c>p</c>.
ALIGN_64 float ps[TEST_CASES];

/// \brief Check function.
void
check()
{
    float e = 1e-3;

    for (int i = 0; i < TEST_CASES; i++)
    {
        float diff_d = abs(ds[i] - ds_orig[i]);
        float diff_u = abs(us[i] - us_orig[i]);
        // Do not check components v and w.
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
void
clean()
{
    for (int i = 0; i < TEST_CASES; i++)
    {
        ds[i] = 0.0;
        us[i] = 0.0;
        vs[i] = 0.0;
        ws_[i] = 0.0;
        ps[i] = 0.0;
    }
}

/// \brief Run riemann solver test and print information.
///
/// \param[in] solver - solver function
/// \param[in] nt - number of threads
/// \param[in] str - description
///
/// \return
/// Execution time.
double
run(void (*solver)(int,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   float *,
                   int),
    int nt,
    string str)
{
    clean();
    double t_start = omp_get_wtime();

    for (int i = 0; i < INNER_REPEATS; i++)
    {
        solver(TEST_CASES,
               dls, uls, vls, wls, pls,
               drs, urs, vrs, wrs, prs,
               ds, us, vs, ws_, ps,
               nt);
    }

    double t_end = omp_get_wtime();
    check();
    double t_len = t_end - t_start;
    cout << setw(15) << str << " ~ " << t_len << " seconds" << endl;

    return t_len;
}

/// \brief Min value in array.
///
/// \param[in] d - array
/// \param[in] c - element count
///
/// \return
/// Min value.
double
array_min(double *d,
          int c)
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
///
/// \parameter[in] argc - arguments count
/// \parameter[in] argv - arguments
int
main(int argc,
     char **argv)
{
    double times[REPEATS];
    double times_opt[REPEATS];
    int dls_len = ARR_LEN(dls);
    int drs_len = ARR_LEN(drs);
    int ds_orig_len = ARR_LEN(ds_orig);
    int uls_len = ARR_LEN(uls);
    int urs_len = ARR_LEN(urs);
    int us_orig_len = ARR_LEN(us_orig);
    int pls_len = ARR_LEN(pls);
    int prs_len = ARR_LEN(prs);
    int ps_orig_len = ARR_LEN(ps_orig);

    // Threads count processing.
    int nt = 1;
    if (argc > 1)
    {
        nt = atoi(argv[1]);
    }

    // Check if data elements count is enough for run.
    if ((TEST_CASES > dls_len)
        || (TEST_CASES > drs_len)
        || (TEST_CASES > ds_orig_len)
        || (TEST_CASES > uls_len)
        || (TEST_CASES > urs_len)
        || (TEST_CASES > us_orig_len)
        || (TEST_CASES > pls_len)
        || (TEST_CASES > prs_len)
        || (TEST_CASES > ps_orig_len))
    {
        cout << "error : not enough data for run : TEST_CASES = " << TEST_CASES
             << ", data = { " << dls_len << ", " << drs_len << ", " << ds_orig_len
             << ", " << uls_len << ", " << urs_len << ", " << us_orig_len
             << ", " << pls_len << ", " << prs_len << ", " << ps_orig_len << " }" << endl;
        exit(1);
    }

#ifdef INTEL

    {
        // Check alignment.
        unsigned long dls_a = (unsigned long)&dls[0];
        unsigned long uls_a = (unsigned long)&uls[0];
        unsigned long vls_a = (unsigned long)&vls[0];
        unsigned long wls_a = (unsigned long)&wls[0];
        unsigned long pls_a = (unsigned long)&pls[0];
        unsigned long drs_a = (unsigned long)&drs[0];
        unsigned long urs_a = (unsigned long)&urs[0];
        unsigned long vrs_a = (unsigned long)&vrs[0];
        unsigned long wrs_a = (unsigned long)&wrs[0];
        unsigned long prs_a = (unsigned long)&prs[0];
        unsigned long ds_orig_a = (unsigned long)&ds_orig[0];
        unsigned long us_orig_a = (unsigned long)&us_orig[0];
        unsigned long ps_orig_a = (unsigned long)&ps_orig[0];
        unsigned long ds_a = (unsigned long)&ds[0];
        unsigned long us_a = (unsigned long)&us[0];
        unsigned long vs_a = (unsigned long)&vs[0];
        unsigned long ws_a = (unsigned long)&ws_[0];
        unsigned long ps_a = (unsigned long)&ps[0];

        if (((dls_a | uls_a | vls_a | wls_a | pls_a
              | drs_a | urs_a | vrs_a | wrs_a | prs_a
              | ds_orig_a | us_orig_a | ps_orig_a
              | ds_a | us_a | vs_a | ws_a | ps_a ) & 0x3F) != 0x0)
        {
            cout << "wrong arrays alignment : " << hex
                 << &dls[0] << ", " << &uls[0] << ", " << &vls[0] << ", " << &wls[0] << ", " << &pls[0] << ", "
                 << &drs[0] << ", " << &urs[0] << ", " << &vrs[0] << ", " << &wrs[0] << ", " << &prs[0] << ", "
                 << &ds_orig[0] << ", " << &us_orig[0] << ", " << &ps_orig[0] << ", "
                 << &ds[0] << ", " << &us[0] << ", " << &vs[0] << ", " << &ws_[0] << ", " << &ps[0]
                 << endl;
            exit(1);
        }
    }

#endif

    cout << "test begin : " << TEST_CASES << " test cases" << endl;

#if !defined(OPENMP_CHUNKS) && !defined(OPENMP_INTERLEAVE) && !defined(OPENMP_RACE) \
    && !defined(OPENMP_CHUNKS_AUTO) && !defined(OPENMP_INTERLEAVE_AUTO) && !defined(OPENMP_RACE_AUTO)
#error "unknown form OpenMP distribution organization"
#endif

#ifdef OPENMP_CHUNKS
    cout << "num_threads = 1 vs " << nt << " (OPENMP_CHUNKS)" << endl;
#endif // OPENMP_CHUNKS

#ifdef OPENMP_INTERLEAVE
    cout << "num_threads = 1 vs " << nt << " (OPENMP_INTERLEAVE)" << endl;
#endif // OPENMP_INTERLEAVE

#ifdef OPENMP_RACE
    cout << "num_threads = 1 vs " << nt << " (OPENMP_RACE)" << endl;
#endif // OPENMP_RACE

#ifdef OPENMP_CHUNKS_AUTO
    cout << "num_threads = 1 vs " << nt << " (OPENMP_CHUNKS_AUTO)" << endl;
#endif // OPENMP_CHUNKS_AUTO

#ifdef OPENMP_INTERLEAVE_AUTO
    cout << "num_threads = 1 vs " << nt << " (OPENMP_INTERLEAVE_AUTO)" << endl;
#endif // OPENMP_INTERLEAVE_AUTO

#ifdef OPENMP_RACE_AUTO
    cout << "num_threads = 1 vs " << nt << " (OPENMP_RACE_AUTO)" << endl;
#endif // OPENMP_RACE_AUTO

    for (int i = 0; i < REPEATS; i++)
    {
        times[i] = run(riemann_n_s, 1, "n_s");
    }

    for (int cur_nt = 1; cur_nt <= nt; cur_nt++)
    {
        cout << "----------" << endl;

        for (int i = 0; i < REPEATS; i++)
        {
            times_opt[i] = run(riemann_n_v, cur_nt, "n_v");
        }

        double min_time = array_min(times, REPEATS);
        double min_time_opt = array_min(times_opt, REPEATS);
        double time_reduce = ((min_time - min_time_opt) / min_time) * 100.0;
        double speedup_x = min_time / min_time_opt;
    
        cout << "test done : "
             << "nt = " << cur_nt
             << ", min_time = " << min_time
             << ", min_time_opt = " << min_time_opt
             << ", time_reduce = " << setprecision(2) << time_reduce
             << "%, speedup_x = " << setprecision(3) << speedup_x << endl;
    }

    return 0;
}
