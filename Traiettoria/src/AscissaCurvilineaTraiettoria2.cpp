//
// File: AscissaCurvilineaTraiettoria2.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "AscissaCurvilineaTraiettoria2.h"
#include "Traj_emxutil.h"
#include "norm.h"
#include "PVATraiettoria2.h"

// Function Definitions

//
// Arguments    : double T
//                const double qi[3]
//                const double qf[3]
//                double tf
//                emxArray_real_T *ps
//                emxArray_real_T *dps
//                emxArray_real_T *ddps
// Return Type  : void
//
void AscissaCurvilineaTraiettoria2(double T, const double qi[3], const double
  qf[3], double tf, emxArray_real_T *ps, emxArray_real_T *dps, emxArray_real_T
  *ddps)
{
  double b_qf[3];
  int i;
  emxArray_real_T *s;
  emxArray_real_T *ds;
  emxArray_real_T *dds;
  double y;
  int tmp;
  emxArray_real_T *b_ps;
  emxArray_real_T *b_dps;
  emxArray_real_T *b_ddps;
  int loop_ub;
  int i1;
  for (i = 0; i < 3; i++) {
    b_qf[i] = qf[i] - qi[i];
  }

  emxInit_real_T(&s, 2);
  emxInit_real_T(&ds, 2);
  emxInit_real_T(&dds, 2);
  PVATraiettoria2(T, norm(b_qf), tf, s, ds, dds);

  //  Posizione
  //  definizione del primo elemento della matrice su cui effettuare la concatenazione 
  for (i = 0; i < 3; i++) {
    b_qf[i] = qf[i] - qi[i];
  }

  y = s->data[0] / norm(b_qf);
  i = ps->size[0] * ps->size[1];
  ps->size[0] = 3;
  ps->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)ps, i, (int)sizeof(double));
  for (i = 0; i < 3; i++) {
    ps->data[i] = qi[i] + y * (qf[i] - qi[i]);
  }

  for (i = 0; i < 3; i++) {
    b_qf[i] = qf[i] - qi[i];
  }

  y = ds->data[0] / norm(b_qf);
  i = dps->size[0] * dps->size[1];
  dps->size[0] = 3;
  dps->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)dps, i, (int)sizeof(double));
  for (i = 0; i < 3; i++) {
    dps->data[i] = y * (qf[i] - qi[i]);
  }

  for (i = 0; i < 3; i++) {
    b_qf[i] = qf[i] - qi[i];
  }

  y = dds->data[0] / norm(b_qf);
  i = ddps->size[0] * ddps->size[1];
  ddps->size[0] = 3;
  ddps->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)ddps, i, (int)sizeof(double));
  for (i = 0; i < 3; i++) {
    ddps->data[i] = y * (qf[i] - qi[i]);
  }

  tmp = 1;
  emxInit_real_T(&b_ps, 2);
  emxInit_real_T(&b_dps, 2);
  emxInit_real_T(&b_ddps, 2);
  while (tmp - 1 <= s->size[1] - 2) {
    for (i = 0; i < 3; i++) {
      b_qf[i] = qf[i] - qi[i];
    }

    y = s->data[s->size[0] * tmp] / norm(b_qf);
    i = b_ps->size[0] * b_ps->size[1];
    b_ps->size[0] = 3;
    b_ps->size[1] = ps->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_ps, i, (int)sizeof(double));
    loop_ub = ps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        b_ps->data[i1 + b_ps->size[0] * i] = ps->data[i1 + ps->size[0] * i];
      }
    }

    for (i = 0; i < 3; i++) {
      b_ps->data[i + b_ps->size[0] * ps->size[1]] = qi[i] + y * (qf[i] - qi[i]);
    }

    i = ps->size[0] * ps->size[1];
    ps->size[0] = 3;
    ps->size[1] = b_ps->size[1];
    emxEnsureCapacity((emxArray__common *)ps, i, (int)sizeof(double));
    loop_ub = b_ps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        ps->data[i1 + ps->size[0] * i] = b_ps->data[i1 + b_ps->size[0] * i];
      }
    }

    for (i = 0; i < 3; i++) {
      b_qf[i] = qf[i] - qi[i];
    }

    y = ds->data[ds->size[0] * tmp] / norm(b_qf);
    i = b_dps->size[0] * b_dps->size[1];
    b_dps->size[0] = 3;
    b_dps->size[1] = dps->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_dps, i, (int)sizeof(double));
    loop_ub = dps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        b_dps->data[i1 + b_dps->size[0] * i] = dps->data[i1 + dps->size[0] * i];
      }
    }

    for (i = 0; i < 3; i++) {
      b_dps->data[i + b_dps->size[0] * dps->size[1]] = y * (qf[i] - qi[i]);
    }

    i = dps->size[0] * dps->size[1];
    dps->size[0] = 3;
    dps->size[1] = b_dps->size[1];
    emxEnsureCapacity((emxArray__common *)dps, i, (int)sizeof(double));
    loop_ub = b_dps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        dps->data[i1 + dps->size[0] * i] = b_dps->data[i1 + b_dps->size[0] * i];
      }
    }

    for (i = 0; i < 3; i++) {
      b_qf[i] = qf[i] - qi[i];
    }

    y = dds->data[dds->size[0] * tmp] / norm(b_qf);
    i = b_ddps->size[0] * b_ddps->size[1];
    b_ddps->size[0] = 3;
    b_ddps->size[1] = ddps->size[1] + 1;
    emxEnsureCapacity((emxArray__common *)b_ddps, i, (int)sizeof(double));
    loop_ub = ddps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        b_ddps->data[i1 + b_ddps->size[0] * i] = ddps->data[i1 + ddps->size[0] *
          i];
      }
    }

    for (i = 0; i < 3; i++) {
      b_ddps->data[i + b_ddps->size[0] * ddps->size[1]] = y * (qf[i] - qi[i]);
    }

    i = ddps->size[0] * ddps->size[1];
    ddps->size[0] = 3;
    ddps->size[1] = b_ddps->size[1];
    emxEnsureCapacity((emxArray__common *)ddps, i, (int)sizeof(double));
    loop_ub = b_ddps->size[1];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < 3; i1++) {
        ddps->data[i1 + ddps->size[0] * i] = b_ddps->data[i1 + b_ddps->size[0] *
          i];
      }
    }

    tmp++;
  }

  emxFree_real_T(&b_ddps);
  emxFree_real_T(&b_dps);
  emxFree_real_T(&b_ps);
  emxFree_real_T(&dds);
  emxFree_real_T(&ds);
  emxFree_real_T(&s);

  // Grafici Posizione  - Velocità - Accelerazione
  //      if qi(1)~=qf(1)
  //          figure(1)
  //          subplot(1,3,1),plot(0:T:tf,ps(1,1:end)), title('Posizione x'), grid on, axis square 
  //          subplot(1,3,2),plot(0:T:tf,dps(1,1:end)), title('Velocità x'), grid on, axis square 
  //          subplot(1,3,3),plot(0:T:tf,ddps(1,1:end)), title('Accelerazione x'), grid on, axis square 
  //      end
  //
  //      if qi(2)~=qf(2)
  //          figure(2)
  //          subplot(1,3,1),plot(0:T:tf,ps(2,1:end)), title('Posizione y'), grid on, axis square 
  //          subplot(1,3,2),plot(0:T:tf,dps(2,1:end)), title('Velocità y'), grid on, axis square 
  //          subplot(1,3,3),plot(0:T:tf,ddps(2,1:end)), title('Accelerazione y'), grid on, axis square 
  //      end
  //
  //      if qi(3)~=qf(3)
  //          figure(3)
  //          subplot(1,3,1),plot(0:T:tf,ps(3,1:end)), title('Posizione z'), grid on, axis square 
  //          subplot(1,3,2),plot(0:T:tf,dps(3,1:end)), title('Velocità z'), grid on, axis square 
  //          subplot(1,3,3),plot(0:T:tf,ddps(3,1:end)), title('Accelerazione z'), grid on, axis square 
  //      end
}

//
// File trailer for AscissaCurvilineaTraiettoria2.cpp
//
// [EOF]
//
