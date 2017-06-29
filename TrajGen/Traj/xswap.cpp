//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xswap.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "xswap.h"

// Function Definitions

//
// Arguments    : int n
//                double x_data[]
//                int ix0
//                int iy0
// Return Type  : void
//
void xswap(int n, double x_data[], int ix0, int iy0)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x_data[ix];
    x_data[ix] = x_data[iy];
    x_data[iy] = temp;
    ix++;
    iy++;
  }
}

//
// File trailer for xswap.cpp
//
// [EOF]
//
