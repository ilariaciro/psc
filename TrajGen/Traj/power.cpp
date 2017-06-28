//
// File: power.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "power.h"
#include "Traj_emxutil.h"

// Function Definitions

//
// Arguments    : const emxArray_real_T *a
//                emxArray_real_T *y
// Return Type  : void
//
void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int k;
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= a->size[1]; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

//
// File trailer for power.cpp
//
// [EOF]
//
