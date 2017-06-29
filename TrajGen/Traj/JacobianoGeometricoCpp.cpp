//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: JacobianoGeometricoCpp.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "JacobianoGeometricoCpp.h"

// Function Definitions

//
// calcolo dello Jacobiano
// Arguments    : const double p[3]
//                const double A[144]
//                double Jg[36]
// Return Type  : void
//
void JacobianoGeometricoCpp(const double p[3], const double A[144], double Jg[36])
{
  int i3;
  int i4;
  double Ab1[16];
  int i5;
  double Ab2[16];
  double Ab3[16];
  double Ab4[16];
  double Ab5[16];
  double b[3];
  double b_b[3];
  double c_b[3];
  double d_b[3];
  double e_b[3];
  double f_b[3];
  double b_A[3];
  double b_Ab1[3];
  double b_Ab2[3];
  double b_Ab3[3];
  double b_Ab4[3];
  double b_Ab5[3];
  for (i3 = 0; i3 < 4; i3++) {
    for (i4 = 0; i4 < 4; i4++) {
      Ab1[i3 + (i4 << 2)] = 0.0;
      for (i5 = 0; i5 < 4; i5++) {
        Ab1[i3 + (i4 << 2)] += A[i3 + (i5 << 2)] * A[i5 + ((4 + i4) << 2)];
      }
    }

    for (i4 = 0; i4 < 4; i4++) {
      Ab2[i3 + (i4 << 2)] = 0.0;
      for (i5 = 0; i5 < 4; i5++) {
        Ab2[i3 + (i4 << 2)] += Ab1[i3 + (i5 << 2)] * A[i5 + ((8 + i4) << 2)];
      }
    }

    for (i4 = 0; i4 < 4; i4++) {
      Ab3[i3 + (i4 << 2)] = 0.0;
      for (i5 = 0; i5 < 4; i5++) {
        Ab3[i3 + (i4 << 2)] += Ab2[i3 + (i5 << 2)] * A[i5 + ((12 + i4) << 2)];
      }
    }

    for (i4 = 0; i4 < 4; i4++) {
      Ab4[i3 + (i4 << 2)] = 0.0;
      for (i5 = 0; i5 < 4; i5++) {
        Ab4[i3 + (i4 << 2)] += Ab3[i3 + (i5 << 2)] * A[i5 + ((16 + i4) << 2)];
      }
    }

    for (i4 = 0; i4 < 4; i4++) {
      Ab5[i3 + (i4 << 2)] = 0.0;
      for (i5 = 0; i5 < 4; i5++) {
        Ab5[i3 + (i4 << 2)] += Ab4[i3 + (i5 << 2)] * A[i5 + ((20 + i4) << 2)];
      }
    }
  }

  //  Vettore z: matrice di rototraslazione che indicano gli assi di rotazione
  //  Vettori p: matrice di rototraslazione che indicano la posizione dei giunti 
  //  che è utile solo quando il giunto è rotoidale
  //  Jacobiano Geometrico
  for (i3 = 0; i3 < 3; i3++) {
    b[i3] = p[i3] - A[12 + i3];
    b_b[i3] = p[i3] - Ab1[12 + i3];
    c_b[i3] = p[i3] - Ab2[12 + i3];
    d_b[i3] = p[i3] - Ab3[12 + i3];
    e_b[i3] = p[i3] - Ab4[12 + i3];
    f_b[i3] = p[i3] - Ab5[12 + i3];
  }

  b_A[0] = A[9] * b[2] - A[10] * b[1];
  b_A[1] = A[10] * b[0] - A[8] * b[2];
  b_A[2] = A[8] * b[1] - A[9] * b[0];
  b_Ab1[0] = Ab1[9] * b_b[2] - Ab1[10] * b_b[1];
  b_Ab1[1] = Ab1[10] * b_b[0] - Ab1[8] * b_b[2];
  b_Ab1[2] = Ab1[8] * b_b[1] - Ab1[9] * b_b[0];
  b_Ab2[0] = Ab2[9] * c_b[2] - Ab2[10] * c_b[1];
  b_Ab2[1] = Ab2[10] * c_b[0] - Ab2[8] * c_b[2];
  b_Ab2[2] = Ab2[8] * c_b[1] - Ab2[9] * c_b[0];
  b_Ab3[0] = Ab3[9] * d_b[2] - Ab3[10] * d_b[1];
  b_Ab3[1] = Ab3[10] * d_b[0] - Ab3[8] * d_b[2];
  b_Ab3[2] = Ab3[8] * d_b[1] - Ab3[9] * d_b[0];
  b_Ab4[0] = Ab4[9] * e_b[2] - Ab4[10] * e_b[1];
  b_Ab4[1] = Ab4[10] * e_b[0] - Ab4[8] * e_b[2];
  b_Ab4[2] = Ab4[8] * e_b[1] - Ab4[9] * e_b[0];
  b_Ab5[0] = Ab5[9] * f_b[2] - Ab5[10] * f_b[1];
  b_Ab5[1] = Ab5[10] * f_b[0] - Ab5[8] * f_b[2];
  b_Ab5[2] = Ab5[8] * f_b[1] - Ab5[9] * f_b[0];
  for (i3 = 0; i3 < 3; i3++) {
    Jg[i3] = b_A[i3];
    Jg[6 + i3] = b_Ab1[i3];
    Jg[12 + i3] = b_Ab2[i3];
    Jg[18 + i3] = b_Ab3[i3];
    Jg[24 + i3] = b_Ab4[i3];
    Jg[30 + i3] = b_Ab5[i3];
    Jg[i3 + 3] = A[8 + i3];
    Jg[i3 + 9] = Ab1[8 + i3];
    Jg[i3 + 15] = Ab2[8 + i3];
    Jg[i3 + 21] = Ab3[8 + i3];
    Jg[i3 + 27] = Ab4[8 + i3];
    Jg[i3 + 33] = Ab5[8 + i3];
  }
}

//
// File trailer for JacobianoGeometricoCpp.cpp
//
// [EOF]
//
