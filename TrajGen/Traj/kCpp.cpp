//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: kCpp.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "kCpp.h"

// Function Declarations
static double rt_atan2d_snf(double u0, double u1);

// Function Definitions

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2((double)b_u0, (double)b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

//
// Procedura Denavit-Hartenberg
// alfa       teta      d           a
// q0
// q1
// q2
// q3
// q4
// q5
// Arguments    : const double q[6]
//                const char str[3]
//                double p[3]
//                double phi_data[]
//                int phi_size[2]
//                double R[9]
//                double A[144]
// Return Type  : void
//
void kCpp(const double q[6], const char str[3], double p[3], double phi_data[],
          int phi_size[2], double R[9], double A[144])
{
  double A01[16];
  double A12[16];
  double A23[16];
  double A34[16];
  double A45[16];
  double A56[16];
  int kstr;
  static const signed char iv1[4] = { 0, 0, 0, 1 };

  double a[16];
  double b_a[16];
  double c_a[16];
  double d_a[16];
  double e_a[16];
  double f_a[16];
  int i1;
  int i2;
  static const double g_a[16] = { -1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 4.4462, 0.5, 1.0 };

  double Abe[16];
  static const signed char b[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1 };

  boolean_T b_bool;
  int exitg2;
  static const char cv3[3] = { 'Z', 'Y', 'X' };

  double cos_teta;
  double teta;
  double dv0[3];
  int exitg1;
  static const char cv4[3] = { 'Z', 'Y', 'Z' };

  double sin_teta;
  double r;
  double dv1[3];

  //  Matrice terna base
  //  Matrice rototraslazione manipolare
  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A01[0] = cos(q[0]);
  A01[4] = -sin(q[0]) * 6.123233995736766E-17;
  A01[8] = -sin(q[0]);
  A01[12] = 0.15 * cos(q[0]);
  A01[1] = sin(q[0]);
  A01[5] = cos(q[0]) * 6.123233995736766E-17;
  A01[9] = -(-cos(q[0]));
  A01[13] = 0.15 * sin(q[0]);
  A01[2] = 0.0;
  A01[6] = -1.0;
  A01[10] = 6.123233995736766E-17;
  A01[14] = 0.45;

  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A12[0] = cos(q[1]);
  A12[4] = -sin(q[1]);
  A12[8] = sin(q[1]) * 0.0;
  A12[12] = 0.59 * cos(q[1]);
  A12[1] = sin(q[1]);
  A12[5] = cos(q[1]);
  A12[9] = -cos(q[1]) * 0.0;
  A12[13] = 0.59 * sin(q[1]);
  A12[2] = 0.0;
  A12[6] = 0.0;
  A12[10] = 1.0;
  A12[14] = 0.0;

  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A23[0] = cos(q[2]);
  A23[4] = -sin(q[2]) * 6.123233995736766E-17;
  A23[8] = -sin(q[2]);
  A23[12] = 0.13 * cos(q[2]);
  A23[1] = sin(q[2]);
  A23[5] = cos(q[2]) * 6.123233995736766E-17;
  A23[9] = -(-cos(q[2]));
  A23[13] = 0.13 * sin(q[2]);
  A23[2] = 0.0;
  A23[6] = -1.0;
  A23[10] = 6.123233995736766E-17;
  A23[14] = 0.0;

  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A34[0] = cos(q[3]);
  A34[4] = -sin(q[3]) * 6.123233995736766E-17;
  A34[8] = sin(q[3]);
  A34[12] = 0.0 * cos(q[3]);
  A34[1] = sin(q[3]);
  A34[5] = cos(q[3]) * 6.123233995736766E-17;
  A34[9] = -cos(q[3]);
  A34[13] = 0.0 * sin(q[3]);
  A34[2] = 0.0;
  A34[6] = 1.0;
  A34[10] = 6.123233995736766E-17;
  A34[14] = 0.6471;

  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A45[0] = cos(q[4]);
  A45[4] = -sin(q[4]) * 6.123233995736766E-17;
  A45[8] = -sin(q[4]);
  A45[12] = 0.0 * cos(q[4]);
  A45[1] = sin(q[4]);
  A45[5] = cos(q[4]) * 6.123233995736766E-17;
  A45[9] = -(-cos(q[4]));
  A45[13] = 0.0 * sin(q[4]);
  A45[2] = 0.0;
  A45[6] = -1.0;
  A45[10] = 6.123233995736766E-17;
  A45[14] = 0.0;

  //  Matrice Denavit - Hartenberg
  //  DH=[alfa       teta     d           a]
  //  Matrice rototraslazione - Cinematica diretta
  //  A = [   cos(teta)   -sin(teta)*cos(alfa)   sin(teta)*sin(alfa)    a*cos(teta); 
  //          sin(teta)   cos(teta)*cos(alfa)    -cos(teta)*sin(alfa)   a*sin(teta); 
  //          0            sin(alfa)               cos(alfa)            d          ; 
  //          0            0                        0                   1
  //      ];
  A56[0] = cos(q[5]);
  A56[4] = -sin(q[5]);
  A56[8] = sin(q[5]) * 0.0;
  A56[12] = 0.0 * cos(q[5]);
  A56[1] = sin(q[5]);
  A56[5] = cos(q[5]);
  A56[9] = -cos(q[5]) * 0.0;
  A56[13] = 0.0 * sin(q[5]);
  A56[2] = 0.0;
  A56[6] = 0.0;
  A56[10] = 1.0;
  A56[14] = 0.095;
  for (kstr = 0; kstr < 4; kstr++) {
    A01[3 + (kstr << 2)] = iv1[kstr];
    A12[3 + (kstr << 2)] = iv1[kstr];
    A23[3 + (kstr << 2)] = iv1[kstr];
    A34[3 + (kstr << 2)] = iv1[kstr];
    A45[3 + (kstr << 2)] = iv1[kstr];
    A56[3 + (kstr << 2)] = iv1[kstr];
  }

  for (kstr = 0; kstr < 4; kstr++) {
    for (i1 = 0; i1 < 4; i1++) {
      a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        a[kstr + (i1 << 2)] += g_a[kstr + (i2 << 2)] * A01[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      b_a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        b_a[kstr + (i1 << 2)] += a[kstr + (i2 << 2)] * A12[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      c_a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        c_a[kstr + (i1 << 2)] += b_a[kstr + (i2 << 2)] * A23[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      d_a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        d_a[kstr + (i1 << 2)] += c_a[kstr + (i2 << 2)] * A34[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      e_a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        e_a[kstr + (i1 << 2)] += d_a[kstr + (i2 << 2)] * A45[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      f_a[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        f_a[kstr + (i1 << 2)] += e_a[kstr + (i2 << 2)] * A56[i2 + (i1 << 2)];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      Abe[kstr + (i1 << 2)] = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        Abe[kstr + (i1 << 2)] += f_a[kstr + (i2 << 2)] * (double)b[i2 + (i1 << 2)];
      }
    }
  }

  //  Posizione terna end-effector (metri)
  //  Matrice di rotazione terna end-effector
  for (kstr = 0; kstr < 3; kstr++) {
    p[kstr] = Abe[12 + kstr];
    for (i1 = 0; i1 < 3; i1++) {
      R[i1 + 3 * kstr] = Abe[i1 + (kstr << 2)];
    }
  }

  //  phi=tform2eul(Abe,str); % FUNZIONE MATLAB
  //  E' stata implementata una funzione che calcola gli angoli di Eulero
  phi_size[0] = 1;
  phi_size[1] = 1;
  phi_data[0] = 0.0;
  b_bool = false;
  kstr = 0;
  do {
    exitg2 = 0;
    if (kstr + 1 < 4) {
      if (str[kstr] != cv3[kstr]) {
        exitg2 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (b_bool) {
    cos_teta = sqrt(Abe[0] * Abe[0] + Abe[1] * Abe[1]);
    if (cos_teta == 0.0) {
      teta = rt_atan2d_snf(-Abe[2], 0.0);
      dv0[0] = 0.0;
      dv0[1] = teta;
      dv0[2] = sin(teta) * rt_atan2d_snf(Abe[4], Abe[5]);
      phi_size[0] = 1;
      phi_size[1] = 3;
      for (kstr = 0; kstr < 3; kstr++) {
        phi_data[kstr] = dv0[kstr];
      }
    } else {
      dv0[0] = rt_atan2d_snf(-Abe[1], -Abe[0]);
      dv0[1] = rt_atan2d_snf(-Abe[2], -cos_teta);
      dv0[2] = rt_atan2d_snf(-Abe[6], -Abe[10]);
      phi_size[0] = 1;
      phi_size[1] = 3;
      for (kstr = 0; kstr < 3; kstr++) {
        phi_data[kstr] = dv0[kstr];
      }
    }
  }

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr + 1 < 4) {
      if (str[kstr] != cv4[kstr]) {
        exitg1 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (b_bool) {
    sin_teta = sqrt(Abe[8] * Abe[8] + Abe[9] * Abe[9]);
    if (sin_teta == 0.0) {
      r = rt_atan2d_snf(0.0, Abe[10]);
      dv1[0] = 0.0;
      dv1[1] = -r;
      dv1[2] = cos(-r) * rt_atan2d_snf(-Abe[4], Abe[5]);
      phi_size[0] = 1;
      phi_size[1] = 3;
      for (kstr = 0; kstr < 3; kstr++) {
        phi_data[kstr] = dv1[kstr];
      }
    } else {
      dv1[0] = rt_atan2d_snf(-Abe[9], -Abe[8]);
      dv1[1] = rt_atan2d_snf(-sin_teta, Abe[10]);
      dv1[2] = rt_atan2d_snf(-Abe[6], Abe[2]);
      phi_size[0] = 1;
      phi_size[1] = 3;
      for (kstr = 0; kstr < 3; kstr++) {
        phi_data[kstr] = dv1[kstr];
      }
    }
  }

  //  Struttura matrici di rototraslazione
  //  A=struct('Ab0',Ab0,'A01',A01,'A12',A12,'A23',A23,'A34',A34,'A45', A45,'A56',A56,'A6e',A6e,'Abe',Abe); 
  for (kstr = 0; kstr < 4; kstr++) {
    for (i1 = 0; i1 < 4; i1++) {
      A[i1 + (kstr << 2)] = g_a[i1 + (kstr << 2)];
      A[i1 + ((kstr + 4) << 2)] = A01[i1 + (kstr << 2)];
      A[i1 + ((kstr + 8) << 2)] = A12[i1 + (kstr << 2)];
      A[i1 + ((kstr + 12) << 2)] = A23[i1 + (kstr << 2)];
      A[i1 + ((kstr + 16) << 2)] = A34[i1 + (kstr << 2)];
      A[i1 + ((kstr + 20) << 2)] = A45[i1 + (kstr << 2)];
      A[i1 + ((kstr + 24) << 2)] = A56[i1 + (kstr << 2)];
      A[i1 + ((kstr + 28) << 2)] = b[i1 + (kstr << 2)];
      A[i1 + ((kstr + 32) << 2)] = Abe[i1 + (kstr << 2)];
    }
  }
}

//
// File trailer for kCpp.cpp
//
// [EOF]
//
