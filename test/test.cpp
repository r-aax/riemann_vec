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
using namespace std;

/// \brief Test data.
double data[][9] =
{

#include "test_data.inc"

};

/// \brief Test.
int main()
{
    int test_cases = sizeof(data) / (sizeof(data[0]));
    double d, u, p;
    double e = 1e-4;

    init_gamas();

    cout << "test begin " << test_cases << " test cases" << endl;

    for (int i = 0; i < test_cases; i++)
    {
        riemann(data[i][0], data[i][1], data[i][2],
                data[i][3], data[i][4], data[i][5],
                d, u, p);

        double diff_d = abs(d - data[i][6]);
        double diff_u = abs(u - data[i][7]);
        double diff_p = abs(p - data[i][8]);

        if (!(diff_d < e) && (diff_u < e) && (diff_p < e))
        {
            cerr << "error : " << endl;
            cerr << "  res : " << d << ", " << u << ", " << p << endl;
            cerr << "right : " << data[i][6] << ", " << data[i][7] << ", " << data[i][8] << endl; 
            cerr << " diff : " << diff_d << ", " << diff_u << ", " << diff_p << endl;
            exit(1);
        }
    }

    cout << "test done" << endl;

    return 0;
}
