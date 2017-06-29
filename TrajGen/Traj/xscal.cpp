//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xscal.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "xscal.h"

// Function Definitions

//
// Arguments    : int n
//                double a
//                double x_data[]
//                int ix0
// Return Type  : void
//
void xscal(int n, double a, double x_data[], int ix0)
{
  int i9;
  int k;
  i9 = (ix0 + n) - 1;
  for (k = ix0; k <= i9; k++) {
    x_data[k - 1] *= a;
  }
}

//
// File trailer for xscal.cpp
//
// [EOF]
//
