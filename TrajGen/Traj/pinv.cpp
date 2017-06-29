//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: pinv.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "pinv.h"
#include "svd1.h"

// Function Definitions

//
// Arguments    : const double A_data[]
//                const int A_size[2]
//                double X_data[]
//                int X_size[2]
// Return Type  : void
//
void eml_pinv(const double A_data[], const int A_size[2], double X_data[], int
              X_size[2])
{
  int m;
  int n;
  int k;
  int i7;
  int V_size[2];
  double V_data[36];
  int S_size[2];
  double S_data[36];
  int U_size[2];
  double U_data[36];
  double tol;
  int r;
  int vcol;
  int br;
  int ic;
  int ar;
  int ib;
  int ia;
  int i8;
  m = A_size[0];
  n = A_size[1];
  X_size[0] = (signed char)A_size[1];
  X_size[1] = (signed char)A_size[0];
  k = (signed char)A_size[1] * (signed char)A_size[0];
  for (i7 = 0; i7 < k; i7++) {
    X_data[i7] = 0.0;
  }

  svd(A_data, A_size, U_data, U_size, S_data, S_size, V_data, V_size);
  tol = (double)A_size[0] * S_data[0] * 2.2204460492503131E-16;
  r = 0;
  k = 0;
  while ((k + 1 <= n) && (S_data[k + S_size[0] * k] > tol)) {
    r++;
    k++;
  }

  if (r > 0) {
    vcol = 0;
    for (br = 0; br + 1 <= r; br++) {
      tol = 1.0 / S_data[br + S_size[0] * br];
      i7 = vcol + n;
      for (k = vcol; k + 1 <= i7; k++) {
        V_data[k] *= tol;
      }

      vcol += n;
    }

    k = A_size[1] * (A_size[0] - 1);
    for (vcol = 0; vcol <= k; vcol += n) {
      i7 = vcol + n;
      for (ic = vcol; ic + 1 <= i7; ic++) {
        X_data[ic] = 0.0;
      }
    }

    br = -1;
    for (vcol = 0; vcol <= k; vcol += n) {
      ar = 0;
      br++;
      i7 = (br + m * (r - 1)) + 1;
      for (ib = br; ib + 1 <= i7; ib += m) {
        if (U_data[ib] != 0.0) {
          ia = ar;
          i8 = vcol + n;
          for (ic = vcol; ic + 1 <= i8; ic++) {
            ia++;
            X_data[ic] += U_data[ib] * V_data[ia - 1];
          }
        }

        ar += n;
      }
    }
  }
}

//
// File trailer for pinv.cpp
//
// [EOF]
//
