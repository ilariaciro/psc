//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xrot.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "xrot.h"

// Function Definitions

//
// Arguments    : int n
//                double x_data[]
//                int ix0
//                int iy0
//                double c
//                double s
// Return Type  : void
//
void xrot(int n, double x_data[], int ix0, int iy0, double c, double s)
{
  int ix;
  int iy;
  int k;
  double temp;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = c * x_data[ix] + s * x_data[iy];
    x_data[iy] = c * x_data[iy] - s * x_data[ix];
    x_data[ix] = temp;
    iy++;
    ix++;
  }
}

//
// File trailer for xrot.cpp
//
// [EOF]
//
