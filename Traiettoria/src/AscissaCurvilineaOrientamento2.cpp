//
// File: AscissaCurvilineaOrientamento2.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "AscissaCurvilineaOrientamento2.h"
#include "Traj_emxutil.h"
#include "norm.h"
#include "PVATraiettoria2.h"

// Function Definitions

//
// T=0.001;
//      phii=[pi; 0; -pi];
//      phif=[pi; 0; -pi];
//      tf=30;
//      dqc=[0; 0; 0.11];
// Arguments    : double T
//                const double phii[3]
//                const double phif[3]
//                double tf
//                emxArray_real_T *phie
//                emxArray_real_T *dphie
//                emxArray_real_T *ddphie
// Return Type  : void
//
void AscissaCurvilineaOrientamento2(double T, const double phii[3], const double
  phif[3], double tf, emxArray_real_T *phie, emxArray_real_T *dphie,
  emxArray_real_T *ddphie)
{
  double b_phif[3];
  int nm1d2;
  double phi_f;
  emxArray_real_T *s;
  int i3;
  int n;
  double anew;
  double apnd;
  double ndbl;
  double cdiff;
  int k;
  int tmp;
  emxArray_real_T *b_phie;
  emxArray_real_T *b_dphie;
  emxArray_real_T *b_ddphie;
  emxArray_real_T *ds;
  emxArray_real_T *dds;
  for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
    b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
  }

  phi_f = norm(b_phif);
  emxInit_real_T(&s, 2);
  if (phi_f == 0.0) {
    i3 = phie->size[0] * phie->size[1];
    phie->size[0] = 3;
    phie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)phie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      phie->data[i3] = phii[i3];
    }

    i3 = dphie->size[0] * dphie->size[1];
    dphie->size[0] = 3;
    dphie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)dphie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      dphie->data[i3] = 0.0;
    }

    i3 = ddphie->size[0] * ddphie->size[1];
    ddphie->size[0] = 3;
    ddphie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)ddphie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      ddphie->data[i3] = 0.0;
    }

    if (rtIsNaN(T) || rtIsNaN(tf)) {
      n = 1;
      anew = rtNaN;
      apnd = tf;
    } else if ((T == 0.0) || ((0.0 < tf) && (T < 0.0)) || ((tf < 0.0) && (T >
                 0.0))) {
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

    i3 = s->size[0] * s->size[1];
    s->size[0] = 1;
    s->size[1] = n;
    emxEnsureCapacity((emxArray__common *)s, i3, (int)sizeof(double));
    if (n > 0) {
      s->data[0] = anew;
      if (n > 1) {
        s->data[n - 1] = apnd;
        i3 = n - 1;
        nm1d2 = i3 / 2;
        for (k = 1; k < nm1d2; k++) {
          ndbl = (double)k * T;
          s->data[k] = anew + ndbl;
          s->data[(n - k) - 1] = apnd - ndbl;
        }

        if (nm1d2 << 1 == n - 1) {
          s->data[nm1d2] = (anew + apnd) / 2.0;
        } else {
          ndbl = (double)nm1d2 * T;
          s->data[nm1d2] = anew + ndbl;
          s->data[nm1d2 + 1] = apnd - ndbl;
        }
      }
    }

    tmp = 0;
    emxInit_real_T(&b_phie, 2);
    emxInit_real_T(&b_dphie, 2);
    emxInit_real_T(&b_ddphie, 2);
    while (tmp <= s->size[1] - 2) {
      i3 = b_phie->size[0] * b_phie->size[1];
      b_phie->size[0] = 3;
      b_phie->size[1] = phie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_phie, i3, (int)sizeof(double));
      nm1d2 = phie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_phie->data[k + b_phie->size[0] * i3] = phie->data[k + phie->size[0] *
            i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_phie->data[i3 + b_phie->size[0] * phie->size[1]] = phii[i3];
      }

      i3 = phie->size[0] * phie->size[1];
      phie->size[0] = 3;
      phie->size[1] = b_phie->size[1];
      emxEnsureCapacity((emxArray__common *)phie, i3, (int)sizeof(double));
      nm1d2 = b_phie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          phie->data[k + phie->size[0] * i3] = b_phie->data[k + b_phie->size[0] *
            i3];
        }
      }

      i3 = b_dphie->size[0] * b_dphie->size[1];
      b_dphie->size[0] = 3;
      b_dphie->size[1] = dphie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_dphie, i3, (int)sizeof(double));
      nm1d2 = dphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_dphie->data[k + b_dphie->size[0] * i3] = dphie->data[k + dphie->
            size[0] * i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_dphie->data[i3 + b_dphie->size[0] * dphie->size[1]] = 0.0;
      }

      i3 = dphie->size[0] * dphie->size[1];
      dphie->size[0] = 3;
      dphie->size[1] = b_dphie->size[1];
      emxEnsureCapacity((emxArray__common *)dphie, i3, (int)sizeof(double));
      nm1d2 = b_dphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          dphie->data[k + dphie->size[0] * i3] = b_dphie->data[k + b_dphie->
            size[0] * i3];
        }
      }

      i3 = b_ddphie->size[0] * b_ddphie->size[1];
      b_ddphie->size[0] = 3;
      b_ddphie->size[1] = ddphie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_ddphie, i3, (int)sizeof(double));
      nm1d2 = ddphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_ddphie->data[k + b_ddphie->size[0] * i3] = ddphie->data[k +
            ddphie->size[0] * i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_ddphie->data[i3 + b_ddphie->size[0] * ddphie->size[1]] = 0.0;
      }

      i3 = ddphie->size[0] * ddphie->size[1];
      ddphie->size[0] = 3;
      ddphie->size[1] = b_ddphie->size[1];
      emxEnsureCapacity((emxArray__common *)ddphie, i3, (int)sizeof(double));
      nm1d2 = b_ddphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          ddphie->data[k + ddphie->size[0] * i3] = b_ddphie->data[k +
            b_ddphie->size[0] * i3];
        }
      }

      tmp++;
    }

    emxFree_real_T(&b_ddphie);
    emxFree_real_T(&b_dphie);
    emxFree_real_T(&b_phie);
  } else {
    emxInit_real_T(&ds, 2);
    emxInit_real_T(&dds, 2);
    PVATraiettoria2(T, phi_f, tf, s, ds, dds);

    //  Orientamento
    for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
      b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
    }

    ndbl = s->data[0] / norm(b_phif);
    i3 = phie->size[0] * phie->size[1];
    phie->size[0] = 3;
    phie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)phie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      phie->data[i3] = phii[i3] + ndbl * (phif[i3] - phii[i3]);
    }

    for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
      b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
    }

    ndbl = ds->data[0] / norm(b_phif);
    i3 = dphie->size[0] * dphie->size[1];
    dphie->size[0] = 3;
    dphie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)dphie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      dphie->data[i3] = ndbl * (phif[i3] - phii[i3]);
    }

    for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
      b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
    }

    ndbl = dds->data[0] / norm(b_phif);
    i3 = ddphie->size[0] * ddphie->size[1];
    ddphie->size[0] = 3;
    ddphie->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)ddphie, i3, (int)sizeof(double));
    for (i3 = 0; i3 < 3; i3++) {
      ddphie->data[i3] = ndbl * (phif[i3] - phii[i3]);
    }

    tmp = 1;
    emxInit_real_T(&b_phie, 2);
    emxInit_real_T(&b_dphie, 2);
    emxInit_real_T(&b_ddphie, 2);
    while (tmp - 1 <= s->size[1] - 2) {
      for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
        b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
      }

      ndbl = s->data[s->size[0] * tmp] / norm(b_phif);
      i3 = b_phie->size[0] * b_phie->size[1];
      b_phie->size[0] = 3;
      b_phie->size[1] = phie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_phie, i3, (int)sizeof(double));
      nm1d2 = phie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_phie->data[k + b_phie->size[0] * i3] = phie->data[k + phie->size[0] *
            i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_phie->data[i3 + b_phie->size[0] * phie->size[1]] = phii[i3] + ndbl *
          (phif[i3] - phii[i3]);
      }

      i3 = phie->size[0] * phie->size[1];
      phie->size[0] = 3;
      phie->size[1] = b_phie->size[1];
      emxEnsureCapacity((emxArray__common *)phie, i3, (int)sizeof(double));
      nm1d2 = b_phie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          phie->data[k + phie->size[0] * i3] = b_phie->data[k + b_phie->size[0] *
            i3];
        }
      }

      for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
        b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
      }

      ndbl = ds->data[ds->size[0] * tmp] / norm(b_phif);
      i3 = b_dphie->size[0] * b_dphie->size[1];
      b_dphie->size[0] = 3;
      b_dphie->size[1] = dphie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_dphie, i3, (int)sizeof(double));
      nm1d2 = dphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_dphie->data[k + b_dphie->size[0] * i3] = dphie->data[k + dphie->
            size[0] * i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_dphie->data[i3 + b_dphie->size[0] * dphie->size[1]] = ndbl * (phif[i3]
          - phii[i3]);
      }

      i3 = dphie->size[0] * dphie->size[1];
      dphie->size[0] = 3;
      dphie->size[1] = b_dphie->size[1];
      emxEnsureCapacity((emxArray__common *)dphie, i3, (int)sizeof(double));
      nm1d2 = b_dphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          dphie->data[k + dphie->size[0] * i3] = b_dphie->data[k + b_dphie->
            size[0] * i3];
        }
      }

      for (nm1d2 = 0; nm1d2 < 3; nm1d2++) {
        b_phif[nm1d2] = phif[nm1d2] - phii[nm1d2];
      }

      ndbl = dds->data[dds->size[0] * tmp] / norm(b_phif);
      i3 = b_ddphie->size[0] * b_ddphie->size[1];
      b_ddphie->size[0] = 3;
      b_ddphie->size[1] = ddphie->size[1] + 1;
      emxEnsureCapacity((emxArray__common *)b_ddphie, i3, (int)sizeof(double));
      nm1d2 = ddphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          b_ddphie->data[k + b_ddphie->size[0] * i3] = ddphie->data[k +
            ddphie->size[0] * i3];
        }
      }

      for (i3 = 0; i3 < 3; i3++) {
        b_ddphie->data[i3 + b_ddphie->size[0] * ddphie->size[1]] = ndbl *
          (phif[i3] - phii[i3]);
      }

      i3 = ddphie->size[0] * ddphie->size[1];
      ddphie->size[0] = 3;
      ddphie->size[1] = b_ddphie->size[1];
      emxEnsureCapacity((emxArray__common *)ddphie, i3, (int)sizeof(double));
      nm1d2 = b_ddphie->size[1];
      for (i3 = 0; i3 < nm1d2; i3++) {
        for (k = 0; k < 3; k++) {
          ddphie->data[k + ddphie->size[0] * i3] = b_ddphie->data[k +
            b_ddphie->size[0] * i3];
        }
      }

      tmp++;
    }

    emxFree_real_T(&b_ddphie);
    emxFree_real_T(&b_dphie);
    emxFree_real_T(&b_phie);
    emxFree_real_T(&dds);
    emxFree_real_T(&ds);

    // Grafici Posizione  - Velocità - Accelerazione angoli Eulero
    //          if phii(1)~=phif(1)
    //              figure(1)
    //              subplot(1,3,1),plot(0:T:tf,phie(1,1:end)), title('Position x'), grid on, axis square 
    //              subplot(1,3,2),plot(0:T:tf,dphie(1,1:end)), title('Velocity x'), grid on, axis square 
    //              subplot(1,3,3),plot(0:T:tf,ddphie(1,1:end)), title('Acceleration x'), grid on, axis square 
    //          end
    //
    //          if phii(2)~=phif(2)
    //              figure(2)
    //              subplot(1,3,1),plot(0:T:tf,phie(2,1:end)), title('Posirion y'), grid on, axis square 
    //              subplot(1,3,2),plot(0:T:tf,dphie(2,1:end)), title('Velocity y'), grid on, axis square 
    //              subplot(1,3,3),plot(0:T:tf,ddphie(2,1:end)), title('Acceleration y'), grid on, axis square 
    //          end
    //
    //          if phii(3)~=phif(3)
    //              figure(3)
    //              subplot(1,3,1),plot(0:T:tf,phie(3,1:end)), title('Position z'), grid on, axis square 
    //              subplot(1,3,2),plot(0:T:tf,dphie(3,1:end)), title('Velocity z'), grid on, axis square 
    //              subplot(1,3,3),plot(0:T:tf,ddphie(3,1:end)), title('Acceleration z'), grid on, axis square 
    //          end
  }

  emxFree_real_T(&s);
}

//
// File trailer for AscissaCurvilineaOrientamento2.cpp
//
// [EOF]
//
