//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: Click1PseudoinversaSamplesCpp.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "pinv.h"
#include "inv.h"
#include "JacobianoGeometricoCpp.h"
#include "kCpp.h"

// Function Declarations
static double rt_roundd_snf(double u);

// Function Definitions

//
// Arguments    : double u
// Return Type  : double
//
static double rt_roundd_snf(double u)
{
  double y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

//
// Funzione utilizzata per la conversione in C++
// Arguments    : const double q[6]
//                const double dq[6]
//                const double xd[6]
//                const double dxd[6]
//                double T
//                const char str[3]
//                double K
//                double qf[6]
//                double dqf_data[]
//                int dqf_size[1]
// Return Type  : void
//
void Click1PseudoinversaSamplesCpp(const double q[6], const double dq[6], const
  double xd[6], const double dxd[6], double T, const char str[3], double K,
  double qf[6], double dqf_data[], int dqf_size[1])
{
  double A[144];
  double unusedU0[9];
  int phi_size[2];
  double phi_data[3];
  double p[3];
  int i0;
  double x_data[6];
  int ar;
  double Jg[36];
  int Ja_size[2];
  double Ja_data[36];
  boolean_T b_bool;
  int exitg6;
  static const char cv0[3] = { 'Z', 'Y', 'Z' };

  int exitg5;
  double y[36];
  int ia;
  static const signed char iv0[18] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0 };

  double b_y[36];
  int exitg4;
  static const char cv1[3] = { 'R', 'P', 'Y' };

  boolean_T guard1 = false;
  boolean_T guard2 = false;
  boolean_T guard3 = false;
  int exitg3;
  static const char cv2[3] = { 'Z', 'Y', 'X' };

  int exitg2;
  int exitg1;
  double r[3];
  int k;
  double b_r;
  int b_Ja_size[2];
  int ib;
  int X_size[2];
  double X_data[36];
  double b_xd[6];
  double b[6];
  static const signed char b_b[36] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  int m;
  int ic;
  int br;

  //  Stato manipolatore
  kCpp(q, str, p, phi_data, phi_size, unusedU0, A);
  for (i0 = 0; i0 < 3; i0++) {
    x_data[i0] = p[i0];
  }

  ar = phi_size[1];
  for (i0 = 0; i0 < ar; i0++) {
    x_data[i0 + 3] = phi_data[phi_size[0] * i0];
  }

  JacobianoGeometricoCpp(*(double (*)[3])&x_data[0], A, Jg);
  Ja_size[0] = 1;
  Ja_size[1] = 1;
  Ja_data[0] = 0.0;
  b_bool = false;
  ar = 0;
  do {
    exitg6 = 0;
    if (ar + 1 < 4) {
      if (str[ar] != cv0[ar]) {
        exitg6 = 1;
      } else {
        ar++;
      }
    } else {
      b_bool = true;
      exitg6 = 1;
    }
  } while (exitg6 == 0);

  if (b_bool && ((fabs(x_data[4]) == 0.0) || (fabs(x_data[4]) ==
        3.1415926535897931))) {
  } else {
    b_bool = false;
    ar = 0;
    do {
      exitg5 = 0;
      if (ar + 1 < 4) {
        if (str[ar] != cv0[ar]) {
          exitg5 = 1;
        } else {
          ar++;
        }
      } else {
        b_bool = true;
        exitg5 = 1;
      }
    } while (exitg5 == 0);

    if (b_bool) {
      //  Jacobiano Analitico
      for (i0 = 0; i0 < 6; i0++) {
        for (ia = 0; ia < 3; ia++) {
          y[ia + 6 * i0] = iv0[ia + 3 * i0];
        }
      }

      for (i0 = 0; i0 < 3; i0++) {
        for (ia = 0; ia < 3; ia++) {
          y[(ia + 6 * i0) + 3] = 0.0;
        }
      }

      y[21] = 0.0;
      y[27] = -sin(x_data[3]);
      y[33] = cos(x_data[3]) * sin(x_data[4]);
      y[22] = 0.0;
      y[28] = cos(x_data[3]);
      y[34] = sin(x_data[3]) * sin(x_data[4]);
      y[23] = 1.0;
      y[29] = 0.0;
      y[35] = cos(x_data[4]);
      invNxN(y, b_y);
      for (i0 = 0; i0 < 6; i0++) {
        for (ia = 0; ia < 6; ia++) {
          y[i0 + 6 * ia] = 0.0;
          for (ar = 0; ar < 6; ar++) {
            y[i0 + 6 * ia] += b_y[i0 + 6 * ar] * Jg[ar + 6 * ia];
          }
        }
      }

      Ja_size[0] = 6;
      Ja_size[1] = 6;
      for (i0 = 0; i0 < 6; i0++) {
        for (ia = 0; ia < 6; ia++) {
          Ja_data[ia + 6 * i0] = y[ia + 6 * i0];
        }
      }
    }
  }

  b_bool = false;
  ar = 0;
  do {
    exitg4 = 0;
    if (ar + 1 < 4) {
      if (str[ar] != cv1[ar]) {
        exitg4 = 1;
      } else {
        ar++;
      }
    } else {
      b_bool = true;
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  guard1 = false;
  guard2 = false;
  guard3 = false;
  if (b_bool) {
    guard3 = true;
  } else {
    b_bool = false;
    ar = 0;
    do {
      exitg3 = 0;
      if (ar + 1 < 4) {
        if (str[ar] != cv2[ar]) {
          exitg3 = 1;
        } else {
          ar++;
        }
      } else {
        b_bool = true;
        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (b_bool) {
      guard3 = true;
    } else {
      guard2 = true;
    }
  }

  if (guard3) {
    if ((fabs(x_data[4]) == 1.5707963267948966) || (fabs(x_data[4]) ==
         4.71238898038469)) {
    } else {
      guard2 = true;
    }
  }

  if (guard2) {
    b_bool = false;
    ar = 0;
    do {
      exitg2 = 0;
      if (ar + 1 < 4) {
        if (str[ar] != cv1[ar]) {
          exitg2 = 1;
        } else {
          ar++;
        }
      } else {
        b_bool = true;
        exitg2 = 1;
      }
    } while (exitg2 == 0);

    if (b_bool) {
      guard1 = true;
    } else {
      b_bool = false;
      ar = 0;
      do {
        exitg1 = 0;
        if (ar + 1 < 4) {
          if (str[ar] != cv2[ar]) {
            exitg1 = 1;
          } else {
            ar++;
          }
        } else {
          b_bool = true;
          exitg1 = 1;
        }
      } while (exitg1 == 0);

      if (b_bool) {
        guard1 = true;
      }
    }
  }

  if (guard1) {
    //  T_phi_RPY=[
    //              cos(teta)cos(psi)         -sin(psi)       0;
    //              cos(teta)sin(psi)         cos(psi)        0;
    //              -sin(teta)                0               1;
    //              ];
    //  Jacobiano Analitico RPY
    for (i0 = 0; i0 < 6; i0++) {
      for (ia = 0; ia < 3; ia++) {
        y[ia + 6 * i0] = iv0[ia + 3 * i0];
      }
    }

    for (i0 = 0; i0 < 3; i0++) {
      for (ia = 0; ia < 3; ia++) {
        y[(ia + 6 * i0) + 3] = 0.0;
      }
    }

    y[21] = cos(x_data[4]) * cos(x_data[3]);
    y[27] = -sin(x_data[3]);
    y[33] = 0.0;
    y[22] = cos(x_data[4]) * sin(x_data[3]);
    y[28] = cos(x_data[3]);
    y[34] = 0.0;
    y[23] = -sin(x_data[4]);
    y[29] = 0.0;
    y[35] = 1.0;
    invNxN(y, b_y);
    for (i0 = 0; i0 < 6; i0++) {
      for (ia = 0; ia < 6; ia++) {
        y[i0 + 6 * ia] = 0.0;
        for (ar = 0; ar < 6; ar++) {
          y[i0 + 6 * ia] += b_y[i0 + 6 * ar] * Jg[ar + 6 * ia];
        }
      }
    }

    Ja_size[0] = 6;
    Ja_size[1] = 6;
    for (i0 = 0; i0 < 6; i0++) {
      for (ia = 0; ia < 6; ia++) {
        Ja_data[ia + 6 * i0] = y[ia + 6 * i0];
      }
    }
  }

  // Elimina le righe nulle dello Jacobiano
  //  toDelete = [];
  //  for i = 1 :size(Ja,1)
  //      if(abs(Ja(i,:)) <10^-8)
  //          toDelete = [toDelete,i];
  //      end
  //  end
  //  Ja(toDelete,:) = [];
  //  Errore
  //      if d<1e-4
  //          d=[0;0;0];
  //      else
  for (k = 0; k < 3; k++) {
    b_r = ((xd[3 + k] - x_data[3 + k]) + 3.1415926535897931) /
      6.2831853071795862;
    if (fabs(b_r - rt_roundd_snf(b_r)) <= 2.2204460492503131E-16 * fabs(b_r)) {
      b_r = 0.0;
    } else {
      b_r = (b_r - floor(b_r)) * 6.2831853071795862;
    }

    r[k] = b_r;
  }

  //      end
  // e_o=xd(4:6,i-1)-x(4:6,i-1);
  //  Inversione cinematica
  if (Ja_size[0] < Ja_size[1]) {
    b_Ja_size[0] = Ja_size[1];
    b_Ja_size[1] = Ja_size[0];
    ar = Ja_size[0];
    for (i0 = 0; i0 < ar; i0++) {
      ib = Ja_size[1];
      for (ia = 0; ia < ib; ia++) {
        y[ia + b_Ja_size[0] * i0] = Ja_data[i0 + Ja_size[0] * ia];
      }
    }

    eml_pinv(y, b_Ja_size, Ja_data, Ja_size);
    X_size[0] = Ja_size[1];
    X_size[1] = Ja_size[0];
    ar = Ja_size[0];
    for (i0 = 0; i0 < ar; i0++) {
      ib = Ja_size[1];
      for (ia = 0; ia < ib; ia++) {
        X_data[ia + X_size[0] * i0] = Ja_data[i0 + Ja_size[0] * ia];
      }
    }
  } else {
    eml_pinv(Ja_data, Ja_size, X_data, X_size);
  }

  for (i0 = 0; i0 < 3; i0++) {
    b_xd[i0] = xd[i0] - x_data[i0];
    b_xd[i0 + 3] = r[i0] - 3.1415926535897931;
  }

  for (i0 = 0; i0 < 6; i0++) {
    b_r = 0.0;
    for (ia = 0; ia < 6; ia++) {
      b_r += K * (double)b_b[i0 + 6 * ia] * b_xd[ia];
    }

    b[i0] = dxd[i0] + b_r;
  }

  if (X_size[1] == 1) {
    dqf_size[0] = X_size[0];
    ar = X_size[0];
    for (i0 = 0; i0 < ar; i0++) {
      dqf_data[i0] = 0.0;
      for (ia = 0; ia < 1; ia++) {
        dqf_data[i0] += X_data[i0] * b[0];
      }
    }
  } else {
    k = X_size[1];
    m = X_size[0];
    ar = (signed char)X_size[0];
    dqf_size[0] = (signed char)X_size[0];
    for (i0 = 0; i0 < ar; i0++) {
      dqf_data[i0] = 0.0;
    }

    ar = 0;
    while (ar <= 0) {
      for (ic = 1; ic <= m; ic++) {
        dqf_data[ic - 1] = 0.0;
      }

      ar = m;
    }

    br = 0;
    ar = 0;
    while (ar <= 0) {
      ar = 0;
      i0 = br + k;
      for (ib = br; ib + 1 <= i0; ib++) {
        if (b[ib] != 0.0) {
          ia = ar;
          for (ic = 0; ic + 1 <= m; ic++) {
            ia++;
            dqf_data[ic] += b[ib] * X_data[ia - 1];
          }
        }

        ar += m;
      }

      br += k;
      ar = m;
    }
  }

  for (ar = 0; ar < 6; ar++) {
    qf[ar] = q[ar] + dq[ar] * T;
  }
}

//
// File trailer for Click1PseudoinversaSamplesCpp.cpp
//
// [EOF]
//
