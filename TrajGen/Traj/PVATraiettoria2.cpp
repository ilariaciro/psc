//
// File: PVATraiettoria2.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "PVATraiettoria2.h"
#include "Traj_emxutil.h"
#include "power.h"

// Function Definitions

//
// Controllo velocità crociera inserita
// Arguments    : double T
//                double qf
//                double tf
//                emxArray_real_T *qd
//                emxArray_real_T *dqd
//                emxArray_real_T *ddqd
// Return Type  : void
//
void PVATraiettoria2(double T, double qf, double tf, emxArray_real_T *qd,
                     emxArray_real_T *dqd, emxArray_real_T *ddqd)
{
  double b_min;
  double dqc;
  double ddqc;
  double tc;
  int n;
  double anew;
  double apnd;
  double ndbl;
  double cdiff;
  emxArray_real_T *t_r;
  int i2;
  int nm1d2;
  int k;
  double kd;
  emxArray_real_T *q_r;
  double absa;
  double absb;
  emxArray_real_T *t_c;
  double y;
  double dq_c;
  emxArray_real_T *t_f;
  emxArray_real_T *b_tf;
  emxArray_real_T *q_f;
  unsigned int uv0[2];
  unsigned int uv1[2];
  unsigned int uv2[2];
  b_min = fabs(qf) / tf;
  dqc = (2.0 * fabs(qf) / tf - b_min) / 2.0 + b_min;

  //  Calcolo accelerazione crociera rispetto a velocità e posizione
  ddqc = dqc * dqc / ((0.0 - qf) + dqc * tf);

  //  Tempo crociera
  tc = ((0.0 - qf) + dqc * tf) / dqc;

  //  Definizione traiettoria
  //  Definizione primo tratto traiettoria
  if (rtIsNaN(T) || rtIsNaN(tc)) {
    n = 1;
    anew = rtNaN;
    apnd = tc;
  } else if ((T == 0.0) || ((0.0 < tc) && (T < 0.0)) || ((tc < 0.0) && (T > 0.0)))
  {
    n = 0;
    anew = 0.0;
    apnd = tc;
  } else if (rtIsInf(tc)) {
    n = 1;
    anew = rtNaN;
    apnd = tc;
  } else if (rtIsInf(T)) {
    n = 1;
    anew = 0.0;
    apnd = tc;
  } else {
    anew = 0.0;
    ndbl = floor(tc / T + 0.5);
    apnd = ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tc;
    } else {
      cdiff = tc - apnd;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(tc)) {
      ndbl++;
      apnd = tc;
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

  emxInit_real_T(&t_r, 2);
  i2 = t_r->size[0] * t_r->size[1];
  t_r->size[0] = 1;
  t_r->size[1] = n;
  emxEnsureCapacity((emxArray__common *)t_r, i2, (int)sizeof(double));
  if (n > 0) {
    t_r->data[0] = anew;
    if (n > 1) {
      t_r->data[n - 1] = apnd;
      i2 = n - 1;
      nm1d2 = i2 / 2;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * T;
        t_r->data[k] = anew + kd;
        t_r->data[(n - k) - 1] = apnd - kd;
      }

      if (nm1d2 << 1 == n - 1) {
        t_r->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * T;
        t_r->data[nm1d2] = anew + kd;
        t_r->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  emxInit_real_T(&q_r, 2);

  //  intervallo temporale
  power(t_r, q_r);
  i2 = q_r->size[0] * q_r->size[1];
  q_r->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)q_r, i2, (int)sizeof(double));
  kd = 0.5 * ddqc;
  nm1d2 = q_r->size[0];
  k = q_r->size[1];
  nm1d2 *= k;
  for (i2 = 0; i2 < nm1d2; i2++) {
    q_r->data[i2] *= kd;
  }

  //  posizione
  //  velocità
  //  accelerazione
  //  Definizione secondo tratto traiettoria
  anew = t_r->data[t_r->size[0] * (t_r->size[1] - 1)] + T;
  kd = tf - tc;
  if (rtIsNaN(anew) || rtIsNaN(T) || rtIsNaN(kd)) {
    n = 1;
    anew = rtNaN;
    apnd = kd;
  } else if ((T == 0.0) || ((anew < kd) && (T < 0.0)) || ((kd < anew) && (T >
               0.0))) {
    n = 0;
    apnd = kd;
  } else if (rtIsInf(anew) || rtIsInf(kd)) {
    n = 1;
    anew = rtNaN;
    apnd = kd;
  } else if (rtIsInf(T)) {
    n = 1;
    apnd = kd;
  } else {
    ndbl = floor((kd - anew) / T + 0.5);
    apnd = anew + ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - kd;
    } else {
      cdiff = kd - apnd;
    }

    absa = fabs(anew);
    absb = fabs(kd);
    if ((absa >= absb) || rtIsNaN(absb)) {
      absb = absa;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
      ndbl++;
      apnd = kd;
    } else if (cdiff > 0.0) {
      apnd = anew + (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&t_c, 2);
  i2 = t_c->size[0] * t_c->size[1];
  t_c->size[0] = 1;
  t_c->size[1] = n;
  emxEnsureCapacity((emxArray__common *)t_c, i2, (int)sizeof(double));
  if (n > 0) {
    t_c->data[0] = anew;
    if (n > 1) {
      t_c->data[n - 1] = apnd;
      i2 = n - 1;
      nm1d2 = i2 / 2;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * T;
        t_c->data[k] = anew + kd;
        t_c->data[(n - k) - 1] = apnd - kd;
      }

      if (nm1d2 << 1 == n - 1) {
        t_c->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * T;
        t_c->data[nm1d2] = anew + kd;
        t_c->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  //  intervallo temporale
  y = tc / 2.0;

  //  posizione
  dq_c = ddqc * tc;

  //  velocità
  //  accelerazione
  //  Definizione terzo tratto traiettoria
  anew = t_c->data[t_c->size[0] * (t_c->size[1] - 1)] + T;
  if (rtIsNaN(anew) || rtIsNaN(T) || rtIsNaN(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if ((T == 0.0) || ((anew < tf) && (T < 0.0)) || ((tf < anew) && (T >
               0.0))) {
    n = 0;
    apnd = tf;
  } else if (rtIsInf(anew) || rtIsInf(tf)) {
    n = 1;
    anew = rtNaN;
    apnd = tf;
  } else if (rtIsInf(T)) {
    n = 1;
    apnd = tf;
  } else {
    ndbl = floor((tf - anew) / T + 0.5);
    apnd = anew + ndbl * T;
    if (T > 0.0) {
      cdiff = apnd - tf;
    } else {
      cdiff = tf - apnd;
    }

    absa = fabs(anew);
    absb = fabs(tf);
    if ((absa >= absb) || rtIsNaN(absb)) {
      absb = absa;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
      ndbl++;
      apnd = tf;
    } else if (cdiff > 0.0) {
      apnd = anew + (ndbl - 1.0) * T;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }
  }

  emxInit_real_T(&t_f, 2);
  i2 = t_f->size[0] * t_f->size[1];
  t_f->size[0] = 1;
  t_f->size[1] = n;
  emxEnsureCapacity((emxArray__common *)t_f, i2, (int)sizeof(double));
  if (n > 0) {
    t_f->data[0] = anew;
    if (n > 1) {
      t_f->data[n - 1] = apnd;
      i2 = n - 1;
      nm1d2 = i2 / 2;
      for (k = 1; k < nm1d2; k++) {
        kd = (double)k * T;
        t_f->data[k] = anew + kd;
        t_f->data[(n - k) - 1] = apnd - kd;
      }

      if (nm1d2 << 1 == n - 1) {
        t_f->data[nm1d2] = (anew + apnd) / 2.0;
      } else {
        kd = (double)nm1d2 * T;
        t_f->data[nm1d2] = anew + kd;
        t_f->data[nm1d2 + 1] = apnd - kd;
      }
    }
  }

  emxInit_real_T(&b_tf, 2);

  //  intervallo temporale
  i2 = b_tf->size[0] * b_tf->size[1];
  b_tf->size[0] = 1;
  b_tf->size[1] = t_f->size[1];
  emxEnsureCapacity((emxArray__common *)b_tf, i2, (int)sizeof(double));
  nm1d2 = t_f->size[0] * t_f->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    b_tf->data[i2] = tf - t_f->data[i2];
  }

  emxInit_real_T(&q_f, 2);
  power(b_tf, q_f);
  i2 = q_f->size[0] * q_f->size[1];
  q_f->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)q_f, i2, (int)sizeof(double));
  kd = 0.5 * ddqc;
  nm1d2 = q_f->size[0];
  k = q_f->size[1];
  nm1d2 *= k;
  emxFree_real_T(&b_tf);
  for (i2 = 0; i2 < nm1d2; i2++) {
    q_f->data[i2] = qf - kd * q_f->data[i2];
  }

  //  posizione
  //  velocità
  //  accelerazione
  //  Traiettoria - Velocità - Accelerazione
  kd = ddqc * tc;
  i2 = qd->size[0] * qd->size[1];
  qd->size[0] = 1;
  qd->size[1] = (q_r->size[1] + t_c->size[1]) + q_f->size[1];
  emxEnsureCapacity((emxArray__common *)qd, i2, (int)sizeof(double));
  nm1d2 = q_r->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    qd->data[qd->size[0] * i2] = q_r->data[q_r->size[0] * i2];
  }

  nm1d2 = t_c->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    qd->data[qd->size[0] * (i2 + q_r->size[1])] = kd * (t_c->data[t_c->size[0] *
      i2] - y);
  }

  nm1d2 = q_f->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    qd->data[qd->size[0] * ((i2 + q_r->size[1]) + t_c->size[1])] = q_f->data
      [q_f->size[0] * i2];
  }

  emxFree_real_T(&q_f);
  emxFree_real_T(&q_r);
  for (i2 = 0; i2 < 2; i2++) {
    uv0[i2] = (unsigned int)t_c->size[i2];
  }

  kd = ddqc * tf;
  i2 = dqd->size[0] * dqd->size[1];
  dqd->size[0] = 1;
  dqd->size[1] = (t_r->size[1] + (int)uv0[1]) + t_f->size[1];
  emxEnsureCapacity((emxArray__common *)dqd, i2, (int)sizeof(double));
  nm1d2 = t_r->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    dqd->data[dqd->size[0] * i2] = ddqc * t_r->data[t_r->size[0] * i2];
  }

  nm1d2 = (int)uv0[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    dqd->data[dqd->size[0] * (i2 + t_r->size[1])] = dq_c;
  }

  nm1d2 = t_f->size[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    dqd->data[dqd->size[0] * ((i2 + t_r->size[1]) + (int)uv0[1])] = kd - ddqc *
      t_f->data[t_f->size[0] * i2];
  }

  for (i2 = 0; i2 < 2; i2++) {
    uv0[i2] = (unsigned int)t_r->size[i2];
  }

  emxFree_real_T(&t_r);
  for (i2 = 0; i2 < 2; i2++) {
    uv1[i2] = (unsigned int)t_c->size[i2];
  }

  emxFree_real_T(&t_c);
  for (i2 = 0; i2 < 2; i2++) {
    uv2[i2] = (unsigned int)t_f->size[i2];
  }

  emxFree_real_T(&t_f);
  k = (int)uv0[1];
  i2 = ddqd->size[0] * ddqd->size[1];
  ddqd->size[0] = 1;
  ddqd->size[1] = ((int)uv0[1] + (int)uv1[1]) + (int)uv2[1];
  emxEnsureCapacity((emxArray__common *)ddqd, i2, (int)sizeof(double));
  nm1d2 = (int)uv0[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    ddqd->data[ddqd->size[0] * i2] = ddqc;
  }

  nm1d2 = (int)uv1[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    ddqd->data[ddqd->size[0] * (i2 + k)] = 0.0;
  }

  nm1d2 = (int)uv2[1];
  for (i2 = 0; i2 < nm1d2; i2++) {
    ddqd->data[ddqd->size[0] * ((i2 + k) + (int)uv1[1])] = -ddqc;
  }
}

//
// File trailer for PVATraiettoria2.cpp
//
// [EOF]
//
