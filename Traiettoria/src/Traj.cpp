//
// File: Traj.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "Traj_emxutil.h"
#include "AscissaCurvilineaOrientamento2.h"
#include "AscissaCurvilineaTraiettoria2.h"

// Function Definitions

//
// Definizione dei tempi
// Arguments    : double T
//                double tf
//                double tStop
//                const double p_i[3]
//                const double p_f[3]
//                const double phi_i[3]
//                const double phi_f[3]
//                emxArray_real_T *xd
//                emxArray_real_T *dxd
//                emxArray_real_T *ddxd
// Return Type  : void
//
void Traj(double T, double tf, double tStop, const double p_i[3], const double
          p_f[3], const double phi_i[3], const double phi_f[3], emxArray_real_T *
          xd, emxArray_real_T *dxd, emxArray_real_T *ddxd)
{
  emxArray_real_T *ps;
  emxArray_real_T *dps;
  emxArray_real_T *ddps;
  emxArray_real_T *phie;
  emxArray_real_T *dphie;
  emxArray_real_T *ddphie;
  int n;
  double anew;
  double apnd;
  double ndbl;
  double cdiff;
  emxArray_real_T *r0;
  int i0;
  int nm1d2;
  int k;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *r3;
  emxArray_real_T *r4;
  emxArray_real_T *r5;
  emxArray_real_T *r6;
  emxArray_real_T *b_ps;
  emxArray_real_T *b_phie;
  emxArray_real_T *b_dps;
  emxArray_real_T *b_dphie;
  emxArray_real_T *b_ddps;
  emxArray_real_T *b_ddphie;
  emxInit_real_T(&ps, 2);
  emxInit_real_T(&dps, 2);
  emxInit_real_T(&ddps, 2);
  emxInit_real_T(&phie, 2);
  emxInit_real_T(&dphie, 2);
  emxInit_real_T(&ddphie, 2);

  //  Definizione Traiettoria - Velocità - Accelazione [Posizione - Orientamento] 
  AscissaCurvilineaTraiettoria2(T, p_i, p_f, tf, ps, dps, ddps);
  AscissaCurvilineaOrientamento2(T, phi_i, phi_f, tf, phie, dphie, ddphie);

  // Campioni aggiuntivi posizione
  if (rtIsNaN(T) || rtIsNaN(tStop)) {
    n = 1;
    anew = rtNaN;
    apnd = tStop;
  } else if ((T == 0.0) || ((0.0 < tStop) && (T < 0.0)) || ((tStop < 0.0) && (T >
    0.0))) {
    n = 0;
    anew = 0.0;
    apnd = tStop;
  } else if (rtIsInf(tStop)) {
    n = 1;
    anew = rtNaN;
    apnd = tStop;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tStop;
  } else {
    anew = 0.0;
    ndbl = floor(tStop / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tStop;
    } else {
      cdiff = tStop - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tStop)) {
      ndbl++;
      apnd = tStop;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r0, 2);
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = 1;
  r0->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r0, i0, (int)sizeof(double));
  if (n > 0) {
    r0->data[0] = anew;
    if (n > 1) {
      r0->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r0->data[k] = anew + ndbl;
        r0->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r0->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r0->data[nm1d2] = anew + ndbl;
        r0->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r1, 2);
  i0 = r1->size[0] * r1->size[1];
  r1->size[0] = 1;
  r1->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r1, i0, (int)sizeof(double));
  if (n > 0) {
    r1->data[0] = anew;
    if (n > 1) {
      r1->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r1->data[k] = anew + ndbl;
        r1->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r1->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r1->data[nm1d2] = anew + ndbl;
        r1->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r2, 2);
  i0 = r2->size[0] * r2->size[1];
  r2->size[0] = 1;
  r2->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r2, i0, (int)sizeof(double));
  if (n > 0) {
    r2->data[0] = anew;
    if (n > 1) {
      r2->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r2->data[k] = anew + ndbl;
        r2->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r2->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r2->data[nm1d2] = anew + ndbl;
        r2->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r3, 2);
  i0 = r3->size[0] * r3->size[1];
  r3->size[0] = 1;
  r3->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r3, i0, (int)sizeof(double));
  if (n > 0) {
    r3->data[0] = anew;
    if (n > 1) {
      r3->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r3->data[k] = anew + ndbl;
        r3->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r3->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r3->data[nm1d2] = anew + ndbl;
        r3->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  // Campioni aggiuntivi orientamento
  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r4, 2);
  i0 = r4->size[0] * r4->size[1];
  r4->size[0] = 1;
  r4->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r4, i0, (int)sizeof(double));
  if (n > 0) {
    r4->data[0] = anew;
    if (n > 1) {
      r4->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r4->data[k] = anew + ndbl;
        r4->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r4->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r4->data[nm1d2] = anew + ndbl;
        r4->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r5, 2);
  i0 = r5->size[0] * r5->size[1];
  r5->size[0] = 1;
  r5->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r5, i0, (int)sizeof(double));
  if (n > 0) {
    r5->data[0] = anew;
    if (n > 1) {
      r5->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r5->data[k] = anew + ndbl;
        r5->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r5->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r5->data[nm1d2] = anew + ndbl;
        r5->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  if (rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tf;
  } else if (rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tf;
  } else {
    anew = 0.0;
    ndbl = floor(tf / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tf)) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&r6, 2);
  i0 = r6->size[0] * r6->size[1];
  r6->size[0] = 1;
  r6->size[1] = n;
  emxEnsureCapacity((emxArray__common *)r6, i0, (int)sizeof(double));
  if (n > 0) {
    r6->data[0] = anew;
    if (n > 1) {
      r6->data[n - 1] = apnd;
      i0 = n - 1;
      nm1d2 = i0 / 2;
      for (k = 1; k < nm1d2; k++) {
        ndbl = (double)k * T;
        r6->data[k] = anew + ndbl;
        r6->data[(n - k) - 1] = apnd - ndbl;
      }

      if (nm1d2 << 1 == n - 1) {
        r6->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        ndbl = (double)nm1d2 * T;
        r6->data[nm1d2] = anew + ndbl;
        r6->data[nm1d2 + 1] = apnd - ndbl;
      }
    }
  }

  emxInit_real_T(&b_ps, 2);

  //  Posizione
  // Orientamento
  //  Traiettoria spazio operativo
  i0 = r1->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_ps->size[0] * b_ps->size[1];
  b_ps->size[0] = 3;
  b_ps->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_ps, n, (int)sizeof(double));
  emxFree_real_T(&r1);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_ps->data[n + b_ps->size[0] * k] = ps->data[n + ps->size[0] * (i0 - 1)];
    }
  }

  emxInit_real_T(&b_phie, 2);
  i0 = r4->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_phie->size[0] * b_phie->size[1];
  b_phie->size[0] = 3;
  b_phie->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_phie, n, (int)sizeof(double));
  emxFree_real_T(&r4);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_phie->data[n + b_phie->size[0] * k] = phie->data[n + phie->size[0] * (i0
        - 1)];
    }
  }

  i0 = xd->size[0] * xd->size[1];
  xd->size[0] = 6;
  xd->size[1] = ps->size[1] + b_ps->size[1];
  emxEnsureCapacity((emxArray__common *)xd, i0, (int)sizeof(double));
  nm1d2 = ps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      xd->data[n + xd->size[0] * i0] = ps->data[n + ps->size[0] * i0];
    }
  }

  nm1d2 = b_ps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      xd->data[n + xd->size[0] * (i0 + ps->size[1])] = b_ps->data[n + b_ps->
        size[0] * i0];
    }
  }

  emxFree_real_T(&b_ps);
  emxFree_real_T(&ps);
  nm1d2 = phie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      xd->data[(n + xd->size[0] * i0) + 3] = phie->data[n + phie->size[0] * i0];
    }
  }

  nm1d2 = b_phie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      xd->data[(n + xd->size[0] * (i0 + phie->size[1])) + 3] = b_phie->data[n +
        b_phie->size[0] * i0];
    }
  }

  emxFree_real_T(&b_phie);
  emxFree_real_T(&phie);
  emxInit_real_T(&b_dps, 2);
  i0 = r2->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_dps->size[0] * b_dps->size[1];
  b_dps->size[0] = 3;
  b_dps->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_dps, n, (int)sizeof(double));
  emxFree_real_T(&r2);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_dps->data[n + b_dps->size[0] * k] = dps->data[n + dps->size[0] * (i0 - 1)];
    }
  }

  emxInit_real_T(&b_dphie, 2);
  i0 = r5->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_dphie->size[0] * b_dphie->size[1];
  b_dphie->size[0] = 3;
  b_dphie->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_dphie, n, (int)sizeof(double));
  emxFree_real_T(&r5);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_dphie->data[n + b_dphie->size[0] * k] = dphie->data[n + dphie->size[0] *
        (i0 - 1)];
    }
  }

  i0 = dxd->size[0] * dxd->size[1];
  dxd->size[0] = 6;
  dxd->size[1] = dps->size[1] + b_dps->size[1];
  emxEnsureCapacity((emxArray__common *)dxd, i0, (int)sizeof(double));
  nm1d2 = dps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      dxd->data[n + dxd->size[0] * i0] = dps->data[n + dps->size[0] * i0];
    }
  }

  nm1d2 = b_dps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      dxd->data[n + dxd->size[0] * (i0 + dps->size[1])] = b_dps->data[n +
        b_dps->size[0] * i0];
    }
  }

  emxFree_real_T(&b_dps);
  emxFree_real_T(&dps);
  nm1d2 = dphie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      dxd->data[(n + dxd->size[0] * i0) + 3] = dphie->data[n + dphie->size[0] *
        i0];
    }
  }

  nm1d2 = b_dphie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      dxd->data[(n + dxd->size[0] * (i0 + dphie->size[1])) + 3] = b_dphie->
        data[n + b_dphie->size[0] * i0];
    }
  }

  emxFree_real_T(&b_dphie);
  emxFree_real_T(&dphie);
  emxInit_real_T(&b_ddps, 2);
  i0 = r3->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_ddps->size[0] * b_ddps->size[1];
  b_ddps->size[0] = 3;
  b_ddps->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_ddps, n, (int)sizeof(double));
  emxFree_real_T(&r3);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_ddps->data[n + b_ddps->size[0] * k] = ddps->data[n + ddps->size[0] * (i0
        - 1)];
    }
  }

  emxInit_real_T(&b_ddphie, 2);
  i0 = r6->size[1];
  nm1d2 = r0->size[1] - 1;
  n = b_ddphie->size[0] * b_ddphie->size[1];
  b_ddphie->size[0] = 3;
  b_ddphie->size[1] = nm1d2;
  emxEnsureCapacity((emxArray__common *)b_ddphie, n, (int)sizeof(double));
  emxFree_real_T(&r6);
  emxFree_real_T(&r0);
  for (n = 0; n < 3; n++) {
    for (k = 0; k < nm1d2; k++) {
      b_ddphie->data[n + b_ddphie->size[0] * k] = ddphie->data[n + ddphie->size
        [0] * (i0 - 1)];
    }
  }

  i0 = ddxd->size[0] * ddxd->size[1];
  ddxd->size[0] = 6;
  ddxd->size[1] = ddps->size[1] + b_ddps->size[1];
  emxEnsureCapacity((emxArray__common *)ddxd, i0, (int)sizeof(double));
  nm1d2 = ddps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      ddxd->data[n + ddxd->size[0] * i0] = ddps->data[n + ddps->size[0] * i0];
    }
  }

  nm1d2 = b_ddps->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      ddxd->data[n + ddxd->size[0] * (i0 + ddps->size[1])] = b_ddps->data[n +
        b_ddps->size[0] * i0];
    }
  }

  emxFree_real_T(&b_ddps);
  emxFree_real_T(&ddps);
  nm1d2 = ddphie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      ddxd->data[(n + ddxd->size[0] * i0) + 3] = ddphie->data[n + ddphie->size[0]
        * i0];
    }
  }

  nm1d2 = b_ddphie->size[1];
  for (i0 = 0; i0 < nm1d2; i0++) {
    for (n = 0; n < 3; n++) {
      ddxd->data[(n + ddxd->size[0] * (i0 + ddphie->size[1])) + 3] =
        b_ddphie->data[n + b_ddphie->size[0] * i0];
    }
  }

  emxFree_real_T(&b_ddphie);
  emxFree_real_T(&ddphie);
}

//
// File trailer for Traj.cpp
//
// [EOF]
//
