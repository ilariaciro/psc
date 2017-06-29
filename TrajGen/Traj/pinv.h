//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: pinv.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//
#ifndef __PINV_H__
#define __PINV_H__

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Click1PseudoinversaSamplesCpp_types.h"

// Function Declarations
extern void eml_pinv(const double A_data[], const int A_size[2], double X_data[],
                     int X_size[2]);

#endif

//
// File trailer for pinv.h
//
// [EOF]
//
