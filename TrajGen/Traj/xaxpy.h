//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: xaxpy.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//
#ifndef __XAXPY_H__
#define __XAXPY_H__

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Click1PseudoinversaSamplesCpp_types.h"

// Function Declarations
extern void b_xaxpy(int n, double a, const double x_data[], int ix0, double
                    y_data[], int iy0);
extern void c_xaxpy(int n, double a, const double x_data[], int ix0, double
                    y_data[], int iy0);
extern void xaxpy(int n, double a, int ix0, double y_data[], int iy0);

#endif

//
// File trailer for xaxpy.h
//
// [EOF]
//
