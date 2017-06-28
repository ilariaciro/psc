//
// File: Traj_types.h
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//
#ifndef __TRAJ_TYPES_H__
#define __TRAJ_TYPES_H__

// Include Files
#include "rtwtypes.h"

// Type Definitions
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray__common

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T
#endif

//
// File trailer for Traj_types.h
//
// [EOF]
//
