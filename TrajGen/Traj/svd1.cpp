//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: svd1.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 29-Jun-2017 11:10:26
//

// Include Files
#include "rt_nonfinite.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "svd1.h"
#include "xaxpy.h"
#include "xdotc.h"
#include "xnrm2.h"
#include "xscal.h"
#include "xrot.h"
#include "xrotg.h"
#include "xswap.h"

// Function Definitions

//
// Arguments    : const double A_data[]
//                const int A_size[2]
//                double U_data[]
//                int U_size[2]
//                double S_data[]
//                int S_size[2]
//                double V_data[]
//                int V_size[2]
// Return Type  : void
//
void svd(const double A_data[], const int A_size[2], double U_data[], int
         U_size[2], double S_data[], int S_size[2], double V_data[], int V_size
         [2])
{
  int m;
  int qs;
  int mm;
  double b_A_data[36];
  int n;
  int p;
  int minnp;
  double s_data[6];
  double e_data[6];
  double work_data[6];
  int Vf_size_idx_0;
  double Vf_data[36];
  int nrt;
  int nct;
  int q;
  int iter;
  int nmq;
  boolean_T apply_transform;
  double ztest0;
  int kase;
  int jj;
  double ztest;
  double snorm;
  boolean_T exitg3;
  boolean_T exitg2;
  double f;
  double b;
  double varargin_1[5];
  double mtmp;
  boolean_T exitg1;
  double sqds;
  m = A_size[0];
  qs = A_size[0] * A_size[1];
  for (mm = 0; mm < qs; mm++) {
    b_A_data[mm] = A_data[mm];
  }

  n = A_size[0];
  p = A_size[1];
  if (A_size[0] <= A_size[1]) {
    minnp = A_size[0];
  } else {
    minnp = A_size[1];
  }

  if (A_size[0] + 1 <= A_size[1]) {
    qs = (signed char)(A_size[0] + 1);
  } else {
    qs = (signed char)A_size[1];
  }

  for (mm = 0; mm < qs; mm++) {
    s_data[mm] = 0.0;
  }

  qs = (signed char)A_size[1];
  for (mm = 0; mm < qs; mm++) {
    e_data[mm] = 0.0;
  }

  qs = (signed char)A_size[0];
  for (mm = 0; mm < qs; mm++) {
    work_data[mm] = 0.0;
  }

  U_size[0] = (signed char)A_size[0];
  U_size[1] = (signed char)minnp;
  qs = (signed char)A_size[0] * (signed char)minnp;
  for (mm = 0; mm < qs; mm++) {
    U_data[mm] = 0.0;
  }

  Vf_size_idx_0 = (signed char)A_size[1];
  qs = (signed char)A_size[1] * (signed char)A_size[1];
  for (mm = 0; mm < qs; mm++) {
    Vf_data[mm] = 0.0;
  }

  if (A_size[1] < 2) {
    qs = 0;
  } else {
    qs = A_size[1] - 2;
  }

  if (qs <= A_size[0]) {
    nrt = qs;
  } else {
    nrt = A_size[0];
  }

  if (A_size[0] < 1) {
    qs = 0;
  } else {
    qs = A_size[0] - 1;
  }

  if (qs <= A_size[1]) {
    nct = qs;
  } else {
    nct = A_size[1];
  }

  if (nct >= nrt) {
    mm = nct;
  } else {
    mm = nrt;
  }

  for (q = 1; q <= mm; q++) {
    iter = q + n * (q - 1);
    nmq = n - q;
    apply_transform = false;
    if (q <= nct) {
      ztest0 = xnrm2(nmq + 1, b_A_data, iter);
      if (ztest0 > 0.0) {
        apply_transform = true;
        if (b_A_data[iter - 1] < 0.0) {
          s_data[q - 1] = -ztest0;
        } else {
          s_data[q - 1] = ztest0;
        }

        if (fabs(s_data[q - 1]) >= 1.0020841800044864E-292) {
          ztest0 = 1.0 / s_data[q - 1];
          kase = iter + nmq;
          for (qs = iter; qs <= kase; qs++) {
            b_A_data[qs - 1] *= ztest0;
          }
        } else {
          kase = iter + nmq;
          for (qs = iter; qs <= kase; qs++) {
            b_A_data[qs - 1] /= s_data[q - 1];
          }
        }

        b_A_data[iter - 1]++;
        s_data[q - 1] = -s_data[q - 1];
      } else {
        s_data[q - 1] = 0.0;
      }
    }

    for (jj = q; jj + 1 <= p; jj++) {
      kase = q + n * jj;
      if (apply_transform) {
        ztest0 = xdotc(nmq + 1, b_A_data, iter, b_A_data, kase);
        xaxpy(nmq + 1, -(ztest0 / b_A_data[(q + m * (q - 1)) - 1]), iter,
              b_A_data, kase);
      }

      e_data[jj] = b_A_data[kase - 1];
    }

    if (q <= nct) {
      for (kase = q - 1; kase + 1 <= n; kase++) {
        U_data[kase + U_size[0] * (q - 1)] = b_A_data[kase + m * (q - 1)];
      }
    }

    if (q <= nrt) {
      iter = p - q;
      ztest0 = b_xnrm2(iter, e_data, q + 1);
      if (ztest0 == 0.0) {
        e_data[q - 1] = 0.0;
      } else {
        if (e_data[q] < 0.0) {
          e_data[q - 1] = -ztest0;
        } else {
          e_data[q - 1] = ztest0;
        }

        ztest0 = e_data[q - 1];
        if (fabs(e_data[q - 1]) >= 1.0020841800044864E-292) {
          ztest0 = 1.0 / e_data[q - 1];
          kase = q + iter;
          for (qs = q; qs + 1 <= kase; qs++) {
            e_data[qs] *= ztest0;
          }
        } else {
          kase = q + iter;
          for (qs = q; qs + 1 <= kase; qs++) {
            e_data[qs] /= ztest0;
          }
        }

        e_data[q]++;
        e_data[q - 1] = -e_data[q - 1];
        if (q + 1 <= n) {
          for (kase = q; kase + 1 <= n; kase++) {
            work_data[kase] = 0.0;
          }

          for (jj = q; jj + 1 <= p; jj++) {
            b_xaxpy(nmq, e_data[jj], b_A_data, (q + n * jj) + 1, work_data, q +
                    1);
          }

          for (jj = q; jj + 1 <= p; jj++) {
            c_xaxpy(nmq, -e_data[jj] / e_data[q], work_data, q + 1, b_A_data, (q
                     + n * jj) + 1);
          }
        }
      }

      for (kase = q; kase + 1 <= p; kase++) {
        Vf_data[kase + Vf_size_idx_0 * (q - 1)] = e_data[kase];
      }
    }
  }

  if (A_size[1] <= A_size[0] + 1) {
    m = A_size[1];
  } else {
    m = A_size[0] + 1;
  }

  if (nct < A_size[1]) {
    s_data[nct] = b_A_data[nct + A_size[0] * nct];
  }

  if (A_size[0] < m) {
    s_data[m - 1] = 0.0;
  }

  if (nrt + 1 < m) {
    e_data[nrt] = b_A_data[nrt + A_size[0] * (m - 1)];
  }

  e_data[m - 1] = 0.0;
  if (nct + 1 <= minnp) {
    for (jj = nct; jj + 1 <= minnp; jj++) {
      for (kase = 1; kase <= n; kase++) {
        U_data[(kase + U_size[0] * jj) - 1] = 0.0;
      }

      U_data[jj + U_size[0] * jj] = 1.0;
    }
  }

  for (q = nct - 1; q + 1 > 0; q--) {
    nmq = n - q;
    iter = q + n * q;
    if (s_data[q] != 0.0) {
      for (jj = q + 1; jj + 1 <= minnp; jj++) {
        kase = (q + n * jj) + 1;
        ztest0 = xdotc(nmq, U_data, iter + 1, U_data, kase);
        xaxpy(nmq, -(ztest0 / U_data[iter]), iter + 1, U_data, kase);
      }

      for (kase = q; kase + 1 <= n; kase++) {
        U_data[kase + U_size[0] * q] = -U_data[kase + U_size[0] * q];
      }

      U_data[iter]++;
      for (kase = 1; kase <= q; kase++) {
        U_data[(kase + U_size[0] * q) - 1] = 0.0;
      }
    } else {
      for (kase = 1; kase <= n; kase++) {
        U_data[(kase + U_size[0] * q) - 1] = 0.0;
      }

      U_data[iter] = 1.0;
    }
  }

  for (q = A_size[1] - 1; q + 1 > 0; q--) {
    if ((q + 1 <= nrt) && (e_data[q] != 0.0)) {
      iter = (p - q) - 1;
      kase = (q + p * q) + 2;
      for (jj = q + 1; jj + 1 <= p; jj++) {
        qs = (q + p * jj) + 2;
        ztest0 = xdotc(iter, Vf_data, kase, Vf_data, qs);
        xaxpy(iter, -(ztest0 / Vf_data[kase - 1]), kase, Vf_data, qs);
      }
    }

    for (kase = 1; kase <= p; kase++) {
      Vf_data[(kase + Vf_size_idx_0 * q) - 1] = 0.0;
    }

    Vf_data[q + Vf_size_idx_0 * q] = 1.0;
  }

  for (q = 0; q + 1 <= m; q++) {
    if (s_data[q] != 0.0) {
      ztest = fabs(s_data[q]);
      ztest0 = s_data[q] / ztest;
      s_data[q] = ztest;
      if (q + 1 < m) {
        e_data[q] /= ztest0;
      }

      if (q + 1 <= n) {
        xscal(n, ztest0, U_data, 1 + n * q);
      }
    }

    if ((q + 1 < m) && (e_data[q] != 0.0)) {
      ztest = fabs(e_data[q]);
      ztest0 = ztest / e_data[q];
      e_data[q] = ztest;
      s_data[q + 1] *= ztest0;
      xscal(p, ztest0, Vf_data, 1 + p * (q + 1));
    }
  }

  mm = m;
  iter = 0;
  snorm = 0.0;
  for (kase = 0; kase + 1 <= m; kase++) {
    ztest0 = fabs(s_data[kase]);
    ztest = fabs(e_data[kase]);
    if ((ztest0 >= ztest) || rtIsNaN(ztest)) {
    } else {
      ztest0 = ztest;
    }

    if ((snorm >= ztest0) || rtIsNaN(ztest0)) {
    } else {
      snorm = ztest0;
    }
  }

  while ((m > 0) && (!(iter >= 75))) {
    q = m - 1;
    exitg3 = false;
    while (!(exitg3 || (q == 0))) {
      ztest0 = fabs(e_data[q - 1]);
      if ((ztest0 <= 2.2204460492503131E-16 * (fabs(s_data[q - 1]) + fabs
            (s_data[q]))) || (ztest0 <= 1.0020841800044864E-292) || ((iter > 20)
           && (ztest0 <= 2.2204460492503131E-16 * snorm))) {
        e_data[q - 1] = 0.0;
        exitg3 = true;
      } else {
        q--;
      }
    }

    if (q == m - 1) {
      kase = 4;
    } else {
      qs = m;
      kase = m;
      exitg2 = false;
      while ((!exitg2) && (kase >= q)) {
        qs = kase;
        if (kase == q) {
          exitg2 = true;
        } else {
          ztest0 = 0.0;
          if (kase < m) {
            ztest0 = fabs(e_data[kase - 1]);
          }

          if (kase > q + 1) {
            ztest0 += fabs(e_data[kase - 2]);
          }

          ztest = fabs(s_data[kase - 1]);
          if ((ztest <= 2.2204460492503131E-16 * ztest0) || (ztest <=
               1.0020841800044864E-292)) {
            s_data[kase - 1] = 0.0;
            exitg2 = true;
          } else {
            kase--;
          }
        }
      }

      if (qs == q) {
        kase = 3;
      } else if (qs == m) {
        kase = 1;
      } else {
        kase = 2;
        q = qs;
      }
    }

    switch (kase) {
     case 1:
      f = e_data[m - 2];
      e_data[m - 2] = 0.0;
      for (qs = m - 2; qs + 1 >= q + 1; qs--) {
        ztest0 = s_data[qs];
        xrotg(&ztest0, &f, &ztest, &b);
        s_data[qs] = ztest0;
        if (qs + 1 > q + 1) {
          f = -b * e_data[qs - 1];
          e_data[qs - 1] *= ztest;
        }

        xrot(p, Vf_data, 1 + p * qs, 1 + p * (m - 1), ztest, b);
      }
      break;

     case 2:
      f = e_data[q - 1];
      e_data[q - 1] = 0.0;
      for (qs = q; qs + 1 <= m; qs++) {
        xrotg(&s_data[qs], &f, &ztest, &b);
        f = -b * e_data[qs];
        e_data[qs] *= ztest;
        xrot(n, U_data, 1 + n * qs, 1 + n * (q - 1), ztest, b);
      }
      break;

     case 3:
      varargin_1[0] = fabs(s_data[m - 1]);
      varargin_1[1] = fabs(s_data[m - 2]);
      varargin_1[2] = fabs(e_data[m - 2]);
      varargin_1[3] = fabs(s_data[q]);
      varargin_1[4] = fabs(e_data[q]);
      qs = 1;
      mtmp = varargin_1[0];
      if (rtIsNaN(varargin_1[0])) {
        kase = 2;
        exitg1 = false;
        while ((!exitg1) && (kase < 6)) {
          qs = kase;
          if (!rtIsNaN(varargin_1[kase - 1])) {
            mtmp = varargin_1[kase - 1];
            exitg1 = true;
          } else {
            kase++;
          }
        }
      }

      if (qs < 5) {
        while (qs + 1 < 6) {
          if (varargin_1[qs] > mtmp) {
            mtmp = varargin_1[qs];
          }

          qs++;
        }
      }

      f = s_data[m - 1] / mtmp;
      ztest0 = s_data[m - 2] / mtmp;
      ztest = e_data[m - 2] / mtmp;
      sqds = s_data[q] / mtmp;
      b = ((ztest0 + f) * (ztest0 - f) + ztest * ztest) / 2.0;
      ztest0 = f * ztest;
      ztest0 *= ztest0;
      if ((b != 0.0) || (ztest0 != 0.0)) {
        ztest = sqrt(b * b + ztest0);
        if (b < 0.0) {
          ztest = -ztest;
        }

        ztest = ztest0 / (b + ztest);
      } else {
        ztest = 0.0;
      }

      f = (sqds + f) * (sqds - f) + ztest;
      ztest0 = sqds * (e_data[q] / mtmp);
      for (qs = q + 1; qs < m; qs++) {
        xrotg(&f, &ztest0, &ztest, &b);
        if (qs > q + 1) {
          e_data[qs - 2] = f;
        }

        f = ztest * s_data[qs - 1] + b * e_data[qs - 1];
        e_data[qs - 1] = ztest * e_data[qs - 1] - b * s_data[qs - 1];
        ztest0 = b * s_data[qs];
        s_data[qs] *= ztest;
        xrot(p, Vf_data, 1 + p * (qs - 1), 1 + p * qs, ztest, b);
        s_data[qs - 1] = f;
        xrotg(&s_data[qs - 1], &ztest0, &ztest, &b);
        f = ztest * e_data[qs - 1] + b * s_data[qs];
        s_data[qs] = -b * e_data[qs - 1] + ztest * s_data[qs];
        ztest0 = b * e_data[qs];
        e_data[qs] *= ztest;
        if (qs < n) {
          xrot(n, U_data, 1 + n * (qs - 1), 1 + n * qs, ztest, b);
        }
      }

      e_data[m - 2] = f;
      iter++;
      break;

     default:
      if (s_data[q] < 0.0) {
        s_data[q] = -s_data[q];
        xscal(p, -1.0, Vf_data, 1 + p * q);
      }

      qs = q + 1;
      while ((q + 1 < mm) && (s_data[q] < s_data[qs])) {
        ztest = s_data[q];
        s_data[q] = s_data[qs];
        s_data[qs] = ztest;
        if (q + 1 < p) {
          xswap(p, Vf_data, 1 + p * q, 1 + p * (q + 1));
        }

        if (q + 1 < n) {
          xswap(n, U_data, 1 + n * q, 1 + n * (q + 1));
        }

        q = qs;
        qs++;
      }

      iter = 0;
      m--;
      break;
    }
  }

  for (qs = 0; qs + 1 <= minnp; qs++) {
    e_data[qs] = s_data[qs];
  }

  V_size[0] = (signed char)A_size[1];
  V_size[1] = (signed char)minnp;
  for (qs = 0; qs + 1 <= minnp; qs++) {
    for (kase = 0; kase + 1 <= p; kase++) {
      V_data[kase + V_size[0] * qs] = Vf_data[kase + Vf_size_idx_0 * qs];
    }
  }

  S_size[0] = minnp;
  S_size[1] = minnp;
  qs = minnp * minnp;
  for (mm = 0; mm < qs; mm++) {
    S_data[mm] = 0.0;
  }

  for (qs = 0; qs < minnp; qs++) {
    S_data[qs + minnp * qs] = e_data[qs];
  }
}

//
// File trailer for svd1.cpp
//
// [EOF]
//
