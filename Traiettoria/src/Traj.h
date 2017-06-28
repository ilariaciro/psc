//
// File: Traj.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//
#ifndef __TRAJ_H__
#define __TRAJ_H__

// Include Files
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "Traj_types.h"

// Function Declarations
extern void Traj(double T, double tf, double tStop, const double p_i[3], const
                 double p_f[3], const double phi_i[3], const double phi_f[3],
                 emxArray_real_T *xd, emxArray_real_T *dxd, emxArray_real_T
                 *ddxd);

#endif

//
// File trailer for Traj.h
//
// [EOF]
//
