//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xdotc.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "xdotc.h"

// Function Definitions

//
// Arguments    : int n
//                const double x_data[]
//                int ix0
//                const double y_data[]
//                int iy0
// Return Type  : double
//
double xdotc(int n, const double x_data[], int ix0, const double y_data[], int
             iy0)
{
  double d;
  int ix;
  int iy;
  int k;
  d = 0.0;
  if (n < 1) {
  } else {
    ix = ix0;
    iy = iy0;
    for (k = 1; k <= n; k++) {
      d += x_data[ix - 1] * y_data[iy - 1];
      ix++;
      iy++;
    }
  }

  return d;
}

//
// File trailer for xdotc.cpp
//
// [EOF]
//
