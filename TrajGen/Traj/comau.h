#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#define RT_PI                          3.14159265358979323846
#define RT_PIF                         3.1415927F
#define RT_LN_10                       2.30258509299404568402
#define RT_LN_10F                      2.3025851F
#define RT_LOG10E                      0.43429448190325182765
#define RT_LOG10EF                     0.43429449F
#define RT_E                           2.7182818284590452354
#define RT_EF                          2.7182817F
#define NumBitsPerChar                 8U
#define MAX_int8_T                     ((int8_T)(127))
#define MIN_int8_T                     ((int8_T)(-128))
#define MAX_uint8_T                    ((uint8_T)(255))
#define MIN_uint8_T                    ((uint8_T)(0))
#define MAX_int16_T                    ((int16_T)(32767))
#define MIN_int16_T                    ((int16_T)(-32768))
#define MAX_uint16_T                   ((uint16_T)(65535))
#define MIN_uint16_T                   ((uint16_T)(0))
#define MAX_int32_T                    ((int32_T)(2147483647))
#define MIN_int32_T                    ((int32_T)(-2147483647-1))
#define MAX_uint32_T                   ((uint32_T)(0xFFFFFFFFU))
#define MIN_uint32_T                   ((uint32_T)(0))
#define MAX_int64_T                    ((int64_T)(9223372036854775807LL))
#define MIN_int64_T                    ((int64_T)(-9223372036854775807LL-1LL))
#define MAX_uint64_T                   ((uint64_T)(0xFFFFFFFFFFFFFFFFULL))
#define MIN_uint64_T                   ((uint64_T)(0ULL))
#define TMW_NAME_LENGTH_MAX            64

typedef signed char int8_T;
typedef unsigned char uint8_T;
typedef short int16_T;
typedef unsigned short uint16_T;
typedef int int32_T;
typedef unsigned int uint32_T;
typedef long long int64_T;
typedef unsigned long long uint64_T;
typedef float real32_T;
typedef double real64_T;
typedef double real_T;
typedef double time_T;
typedef unsigned char boolean_T;
typedef int int_T;
typedef unsigned int uint_T;
typedef unsigned long ulong_T;
typedef unsigned long long ulonglong_T;
typedef char char_T;
typedef char_T byte_T;

real_T rtInf;
real_T rtMinusInf;
real_T rtNaN;
real32_T rtInfF;
real32_T rtMinusInfF;
real32_T rtNaNF;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
#include <float.h>
#endif

typedef struct {
	struct {
		uint32_T wordH;
		uint32_T wordL;
	} words;
} BigEndianIEEEDouble;

typedef struct {
	struct {
		uint32_T wordL;
		uint32_T wordH;
	} words;
} LittleEndianIEEEDouble;

typedef struct {
	union {
		real32_T wordLreal;
		uint32_T wordLuint;
	} wordL;
} IEEESingle;

#define CREAL_T

typedef struct {
	real32_T re;
	real32_T im;
} creal32_T;

typedef struct {
	real64_T re;
	real64_T im;
} creal64_T;

typedef struct {
	real_T re;
	real_T im;
} creal_T;

typedef struct {
	int8_T re;
	int8_T im;
} cint8_T;

typedef struct {
	uint8_T re;
	uint8_T im;
} cuint8_T;

typedef struct {
	int16_T re;
	int16_T im;
} cint16_T;

typedef struct {
	uint16_T re;
	uint16_T im;
} cuint16_T;

typedef struct {
	int32_T re;
	int32_T im;
} cint32_T;

typedef struct {
	uint32_T re;
	uint32_T im;
} cuint32_T;

typedef struct {
	int64_T re;
	int64_T im;
} cint64_T;

typedef struct {
	uint64_T re;
	uint64_T im;
} cuint64_T;


#if !defined(__cplusplus) && !defined(__true_false_are_keywords)
#  ifndef false
#   define false                       (0U)
#  endif

#  ifndef true
#   define true                        (1U)
#  endif
#endif

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
#endif

// PROTOTIPI METODI
void b_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0);
double b_xnrm2(int n, const double x_data[], int ix0);
void c_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0);
void eml_pinv(const double A_data[], const int A_size[2], double X_data[], int
		X_size[2]);
void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize);
void emxFree_real_T(emxArray_real_T **pEmxArray);
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
void invNxN(const double x[36], double y[36]);
void JacobianoGeometricoCpp(const double p[3], const double A[144], double Jg[36]);
void kCpp(const double q[6], const char str[3], double p[3], double phi_data[],
		int phi_size[2], double R[9], double A[144]);
double norm(const double x[3]);
void PVATraiettoria2(double T, double qf, double tf, emxArray_real_T *qd,
		emxArray_real_T *dqd, emxArray_real_T *ddqd);
void rt_InitInfAndNaN(size_t realSize);
real_T rtGetInf(void);
real32_T rtGetInfF(void);
real_T rtGetMinusInf(void);
real32_T rtGetMinusInfF(void);
real_T rtGetNaN(void);
real32_T rtGetNaNF(void);
boolean_T rtIsInf(real_T value);
boolean_T rtIsNaN(real_T value);
void svd(const double A_data[], const int A_size[2], double U_data[], int
		U_size[2], double S_data[], int S_size[2], double V_data[], int V_size
		[2]);
void xaxpy(int n, double a, int ix0, double y_data[], int iy0);
double xdotc(int n, const double x_data[], int ix0, const double y_data[], int iy0);
double xnrm2(int n, const double x_data[], int ix0);
void xrot(int n, double x_data[], int ix0, int iy0, double c, double s);
void xrotg(double *a, double *b, double *c, double *s);
void xscal(int n, double a, double x_data[], int ix0);
void xswap(int n, double x_data[], int ix0, int iy0);

// IMPLEMENTAZIONE METODI
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

		// Grafici Posizione  - Velocit� - Accelerazione angoli Eulero
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

}

void Click1PseudoinversaSamplesCpp_initialize()
{
	rt_InitInfAndNaN(8U);
}

void Click1PseudoinversaSamplesCpp_terminate()
{
	// (no terminate code required)
}

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

void invNxN(const double x[36], double y[36])
{
	double A[36];
	int i6;
	signed char ipiv[6];
	int j;
	int c;
	int jBcol;
	int ix;
	double smax;
	int k;
	double s;
	int i;
	int kAcol;
	signed char p[6];
	for (i6 = 0; i6 < 36; i6++) {
		y[i6] = 0.0;
		A[i6] = x[i6];
	}

	for (i6 = 0; i6 < 6; i6++) {
		ipiv[i6] = (signed char)(1 + i6);
	}

	for (j = 0; j < 5; j++) {
		c = j * 7;
		jBcol = 0;
		ix = c;
		smax = fabs(A[c]);
		for (k = 2; k <= 6 - j; k++) {
			ix++;
			s = fabs(A[ix]);
			if (s > smax) {
				jBcol = k - 1;
				smax = s;
			}
		}

		if (A[c + jBcol] != 0.0) {
			if (jBcol != 0) {
				ipiv[j] = (signed char)((j + jBcol) + 1);
				ix = j;
				jBcol += j;
				for (k = 0; k < 6; k++) {
					smax = A[ix];
					A[ix] = A[jBcol];
					A[jBcol] = smax;
					ix += 6;
					jBcol += 6;
				}
			}

			i6 = (c - j) + 6;
			for (i = c + 1; i + 1 <= i6; i++) {
				A[i] /= A[c];
			}
		}

		jBcol = c;
		kAcol = c + 6;
		for (i = 1; i <= 5 - j; i++) {
			smax = A[kAcol];
			if (A[kAcol] != 0.0) {
				ix = c + 1;
				i6 = (jBcol - j) + 12;
				for (k = 7 + jBcol; k + 1 <= i6; k++) {
					A[k] += A[ix] * -smax;
					ix++;
				}
			}

			kAcol += 6;
			jBcol += 6;
		}
	}

	for (i6 = 0; i6 < 6; i6++) {
		p[i6] = (signed char)(1 + i6);
	}

	for (k = 0; k < 5; k++) {
		if (ipiv[k] > 1 + k) {
			jBcol = p[ipiv[k] - 1];
			p[ipiv[k] - 1] = p[k];
			p[k] = (signed char)jBcol;
		}
	}

	for (k = 0; k < 6; k++) {
		c = p[k] - 1;
		y[k + 6 * (p[k] - 1)] = 1.0;
		for (j = k; j + 1 < 7; j++) {
			if (y[j + 6 * c] != 0.0) {
				for (i = j + 1; i + 1 < 7; i++) {
					y[i + 6 * c] -= y[j + 6 * c] * A[i + 6 * j];
				}
			}
		}
	}

	for (j = 0; j < 6; j++) {
		jBcol = 6 * j;
		for (k = 5; k >= 0; k += -1) {
			kAcol = 6 * k;
			if (y[k + jBcol] != 0.0) {
				y[k + jBcol] /= A[k + kAcol];
				for (i = 0; i + 1 <= k; i++) {
					y[i + jBcol] -= y[k + jBcol] * A[i + kAcol];
				}
			}
		}
	}
}

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
	//  che � utile solo quando il giunto � rotoidale
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

	// Angoli di Eulero
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

double norm(const double x[3])
{
	double y;
	double scale;
	int k;
	double absxk;
	double t;
	y = 0.0;
	scale = 2.2250738585072014E-308;
	for (k = 0; k < 3; k++) {
		absxk = fabs(x[k]);
		if (absxk > scale) {
			t = scale / absxk;
			y = 1.0 + y * t * t;
			scale = absxk;
		} else {
			t = absxk / scale;
			y += t * t;
		}
	}

	return scale * sqrt(y);
}

void eml_pinv(const double A_data[], const int A_size[2], double X_data[], int
		X_size[2])
{
	int m;
	int n;
	int k;
	int i7;
	int V_size[2];
	double V_data[36];
	int S_size[2];
	double S_data[36];
	int U_size[2];
	double U_data[36];
	double tol;
	int r;
	int vcol;
	int br;
	int ic;
	int ar;
	int ib;
	int ia;
	int i8;
	m = A_size[0];
	n = A_size[1];
	X_size[0] = (signed char)A_size[1];
	X_size[1] = (signed char)A_size[0];
	k = (signed char)A_size[1] * (signed char)A_size[0];
	for (i7 = 0; i7 < k; i7++) {
		X_data[i7] = 0.0;
	}

	svd(A_data, A_size, U_data, U_size, S_data, S_size, V_data, V_size);
	tol = (double)A_size[0] * S_data[0] * 2.2204460492503131E-16;
	r = 0;
	k = 0;
	while ((k + 1 <= n) && (S_data[k + S_size[0] * k] > tol)) {
		r++;
		k++;
	}

	if (r > 0) {
		vcol = 0;
		for (br = 0; br + 1 <= r; br++) {
			tol = 1.0 / S_data[br + S_size[0] * br];
			i7 = vcol + n;
			for (k = vcol; k + 1 <= i7; k++) {
				V_data[k] *= tol;
			}

			vcol += n;
		}

		k = A_size[1] * (A_size[0] - 1);
		for (vcol = 0; vcol <= k; vcol += n) {
			i7 = vcol + n;
			for (ic = vcol; ic + 1 <= i7; ic++) {
				X_data[ic] = 0.0;
			}
		}

		br = -1;
		for (vcol = 0; vcol <= k; vcol += n) {
			ar = 0;
			br++;
			i7 = (br + m * (r - 1)) + 1;
			for (ib = br; ib + 1 <= i7; ib += m) {
				if (U_data[ib] != 0.0) {
					ia = ar;
					i8 = vcol + n;
					for (ic = vcol; ic + 1 <= i8; ic++) {
						ia++;
						X_data[ic] += U_data[ib] * V_data[ia - 1];
					}
				}

				ar += n;
			}
		}
	}
}

void power(const emxArray_real_T *a, emxArray_real_T *y)
{
	int k;
	k = y->size[0] * y->size[1];
	y->size[0] = 1;
	y->size[1] = a->size[1];
	emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
	for (k = 0; k + 1 <= a->size[1]; k++) {
		y->data[k] = a->data[k] * a->data[k];
	}
}

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

	//  Calcolo accelerazione crociera rispetto a velocit� e posizione
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
	//  velocit�
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

	//  velocit�
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
	//  velocit�
	//  accelerazione
	//  Traiettoria - Velocit� - Accelerazione
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

void rt_InitInfAndNaN(size_t realSize)
{
	(void) (realSize);
	rtNaN = rtGetNaN();
	rtNaNF = rtGetNaNF();
	rtInf = rtGetInf();
	rtInfF = rtGetInfF();
	rtMinusInf = rtGetMinusInf();
	rtMinusInfF = rtGetMinusInfF();
}

boolean_T rtIsInf(real_T value)
{
	return ((value==rtInf || value==rtMinusInf) ? 1U : 0U);
}

boolean_T rtIsInfF(real32_T value)
{
	return(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
}

boolean_T rtIsNaN(real_T value)
{

#if defined(_MSC_VER) && (_MSC_VER <= 1200)

	return _isnan(value)? TRUE:FALSE;

#else

	return (value!=value)? 1U:0U;

#endif

}

boolean_T rtIsNaNF(real32_T value)
{

#if defined(_MSC_VER) && (_MSC_VER <= 1200)

	return _isnan((real_T)value)? true:false;

#else

	return (value!=value)? 1U:0U;

#endif

}

real_T rtGetInf(void)
{
	size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
	real_T inf = 0.0;
	if (bitsPerReal == 32U) {
		inf = rtGetInfF();
	} else {
		uint16_T one = 1U;
		enum {
			LittleEndian,
			BigEndian
		} machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
		switch (machByteOrder) {
		case LittleEndian:
		{
			union {
				LittleEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0x7FF00000U;
			tmpVal.bitVal.words.wordL = 0x00000000U;
			inf = tmpVal.fltVal;
			break;
		}

		case BigEndian:
		{
			union {
				BigEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0x7FF00000U;
			tmpVal.bitVal.words.wordL = 0x00000000U;
			inf = tmpVal.fltVal;
			break;
		}
		}
	}

	return inf;
}

real32_T rtGetInfF(void)
{
	IEEESingle infF;
	infF.wordL.wordLuint = 0x7F800000U;
	return infF.wordL.wordLreal;
}

real_T rtGetMinusInf(void)
{
	size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
	real_T minf = 0.0;
	if (bitsPerReal == 32U) {
		minf = rtGetMinusInfF();
	} else {
		uint16_T one = 1U;
		enum {
			LittleEndian,
			BigEndian
		} machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
		switch (machByteOrder) {
		case LittleEndian:
		{
			union {
				LittleEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0xFFF00000U;
			tmpVal.bitVal.words.wordL = 0x00000000U;
			minf = tmpVal.fltVal;
			break;
		}

		case BigEndian:
		{
			union {
				BigEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0xFFF00000U;
			tmpVal.bitVal.words.wordL = 0x00000000U;
			minf = tmpVal.fltVal;
			break;
		}
		}
	}

	return minf;
}

real32_T rtGetMinusInfF(void)
{
	IEEESingle minfF;
	minfF.wordL.wordLuint = 0xFF800000U;
	return minfF.wordL.wordLreal;
}

real_T rtGetNaN(void)
{
	size_t bitsPerReal = sizeof(real_T) * (NumBitsPerChar);
	real_T nan = 0.0;
	if (bitsPerReal == 32U) {
		nan = rtGetNaNF();
	} else {
		uint16_T one = 1U;
		enum {
			LittleEndian,
			BigEndian
		} machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
		switch (machByteOrder) {
		case LittleEndian:
		{
			union {
				LittleEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0xFFF80000U;
			tmpVal.bitVal.words.wordL = 0x00000000U;
			nan = tmpVal.fltVal;
			break;
		}

		case BigEndian:
		{
			union {
				BigEndianIEEEDouble bitVal;
				real_T fltVal;
			} tmpVal;

			tmpVal.bitVal.words.wordH = 0x7FFFFFFFU;
			tmpVal.bitVal.words.wordL = 0xFFFFFFFFU;
			nan = tmpVal.fltVal;
			break;
		}
		}
	}

	return nan;
}

real32_T rtGetNaNF(void)
{
	IEEESingle nanF = { { 0 } };

	uint16_T one = 1U;
	enum {
		LittleEndian,
		BigEndian
	} machByteOrder = (*((uint8_T *) &one) == 1U) ? LittleEndian : BigEndian;
	switch (machByteOrder) {
	case LittleEndian:
	{
		nanF.wordL.wordLuint = 0xFFC00000U;
		break;
	}

	case BigEndian:
	{
		nanF.wordL.wordLuint = 0x7FFFFFFFU;
		break;
	}
	}

	return nanF.wordL.wordLreal;
}

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
						b_xaxpy(nmq, e_data[jj], b_A_data, (q + n * jj) + 1, work_data, q +1);
					}

					for (jj = q; jj + 1 <= p; jj++) {
						c_xaxpy(nmq, -e_data[jj] / e_data[q], work_data, q + 1, b_A_data, (q + n * jj) + 1);
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

emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size)
{
	emxArray_real_T *emx;
	int numEl;
	int i;
	emxInit_real_T(&emx, numDimensions);
	numEl = 1;
	for (i = 0; i < numDimensions; i++) {
		numEl *= size[i];
		emx->size[i] = size[i];
	}

	emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
	emx->numDimensions = numDimensions;
	emx->allocatedSize = numEl;
	return emx;
}

emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
		size)
{
	emxArray_real_T *emx;
	int numEl;
	int i;
	emxInit_real_T(&emx, numDimensions);
	numEl = 1;
	for (i = 0; i < numDimensions; i++) {
		numEl *= size[i];
		emx->size[i] = size[i];
	}

	emx->data = data;
	emx->numDimensions = numDimensions;
	emx->allocatedSize = numEl;
	emx->canFreeData = false;
	return emx;
}

emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols)
{
	emxArray_real_T *emx;
	int size[2];
	int numEl;
	int i;
	size[0] = rows;
	size[1] = cols;
	emxInit_real_T(&emx, 2);
	numEl = 1;
	for (i = 0; i < 2; i++) {
		numEl *= size[i];
		emx->size[i] = size[i];
	}

	emx->data = data;
	emx->numDimensions = 2;
	emx->allocatedSize = numEl;
	emx->canFreeData = false;
	return emx;
}

emxArray_real_T *emxCreate_real_T(int rows, int cols)
{
	emxArray_real_T *emx;
	int size[2];
	int numEl;
	int i;
	size[0] = rows;
	size[1] = cols;
	emxInit_real_T(&emx, 2);
	numEl = 1;
	for (i = 0; i < 2; i++) {
		numEl *= size[i];
		emx->size[i] = size[i];
	}

	emx->data = (double *)calloc((unsigned int)numEl, sizeof(double));
	emx->numDimensions = 2;
	emx->allocatedSize = numEl;
	return emx;
}

void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
	emxFree_real_T(&emxArray);
}

void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
	emxInit_real_T(pEmxArray, numDimensions);
}

void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize)
{
	int newNumel;
	int i;
	void *newData;
	newNumel = 1;
	for (i = 0; i < emxArray->numDimensions; i++) {
		newNumel *= emxArray->size[i];
	}

	if (newNumel > emxArray->allocatedSize) {
		i = emxArray->allocatedSize;
		if (i < 16) {
			i = 16;
		}

		while (i < newNumel) {
			i <<= 1;
		}

		newData = calloc((unsigned int)i, (unsigned int)elementSize);
		if (emxArray->data != NULL) {
			memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
			if (emxArray->canFreeData) {
				free(emxArray->data);
			}
		}

		emxArray->data = newData;
		emxArray->allocatedSize = i;
		emxArray->canFreeData = true;
	}
}

void emxFree_real_T(emxArray_real_T **pEmxArray)
{
	if (*pEmxArray != (emxArray_real_T *)NULL) {
		if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
			free((void *)(*pEmxArray)->data);
		}

		free((void *)(*pEmxArray)->size);
		free((void *)*pEmxArray);
		*pEmxArray = (emxArray_real_T *)NULL;
	}
}

void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
	emxArray_real_T *emxArray;
	int i;
	*pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
	emxArray = *pEmxArray;
	emxArray->data = (double *)NULL;
	emxArray->numDimensions = numDimensions;
	emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
	emxArray->allocatedSize = 0;
	emxArray->canFreeData = true;
	for (i = 0; i < numDimensions; i++) {
		emxArray->size[i] = 0;
	}
}

void Traj_initialize()
{
	rt_InitInfAndNaN(8U);
}

void Traj_terminate()
{
	// (no terminate code required)
}

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

	//  Definizione Traiettoria - Velocit� - Accelazione [Posizione - Orientamento]
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

void b_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0)
{
	int ix;
	int iy;
	int k;
	if ((n < 1) || (a == 0.0)) {
	} else {
		ix = ix0 - 1;
		iy = iy0 - 1;
		for (k = 0; k < n; k++) {
			y_data[iy] += a * x_data[ix];
			ix++;
			iy++;
		}
	}
}

void c_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0)
{
	int ix;
	int iy;
	int k;
	if ((n < 1) || (a == 0.0)) {
	} else {
		ix = ix0 - 1;
		iy = iy0 - 1;
		for (k = 0; k < n; k++) {
			y_data[iy] += a * x_data[ix];
			ix++;
			iy++;
		}
	}
}

void xaxpy(int n, double a, int ix0, double y_data[], int iy0)
{
	int ix;
	int iy;
	int k;
	if ((n < 1) || (a == 0.0)) {
	} else {
		ix = ix0 - 1;
		iy = iy0 - 1;
		for (k = 0; k < n; k++) {
			y_data[iy] += a * y_data[ix];
			ix++;
			iy++;
		}
	}
}

double xdotc(int n, const double x_data[], int ix0, const double y_data[], int iy0)
{
	double d;
	int ix;
	int iy;
	int k;
	d = 0.0;
	if (n < 1) {
	} else {
		ix = ix0;
		iy = iy0;
		for (k = 1; k <= n; k++) {
			d += x_data[ix - 1] * y_data[iy - 1];
			ix++;
			iy++;
		}
	}

	return d;
}

double b_xnrm2(int n, const double x_data[], int ix0)
{
	double y;
	double scale;
	int kend;
	int k;
	double absxk;
	double t;
	y = 0.0;
	if (n < 1) {
	} else if (n == 1) {
		y = fabs(x_data[ix0 - 1]);
	} else {
		scale = 2.2250738585072014E-308;
		kend = (ix0 + n) - 1;
		for (k = ix0; k <= kend; k++) {
			absxk = fabs(x_data[k - 1]);
			if (absxk > scale) {
				t = scale / absxk;
				y = 1.0 + y * t * t;
				scale = absxk;
			} else {
				t = absxk / scale;
				y += t * t;
			}
		}

		y = scale * sqrt(y);
	}

	return y;
}

double xnrm2(int n, const double x_data[], int ix0)
{
	double y;
	double scale;
	int kend;
	int k;
	double absxk;
	double t;
	y = 0.0;
	if (n < 1) {
	} else if (n == 1) {
		y = fabs(x_data[ix0 - 1]);
	} else {
		scale = 2.2250738585072014E-308;
		kend = (ix0 + n) - 1;
		for (k = ix0; k <= kend; k++) {
			absxk = fabs(x_data[k - 1]);
			if (absxk > scale) {
				t = scale / absxk;
				y = 1.0 + y * t * t;
				scale = absxk;
			} else {
				t = absxk / scale;
				y += t * t;
			}
		}

		y = scale * sqrt(y);
	}

	return y;
}

void xrot(int n, double x_data[], int ix0, int iy0, double c, double s)
{
	int ix;
	int iy;
	int k;
	double temp;
	ix = ix0 - 1;
	iy = iy0 - 1;
	for (k = 1; k <= n; k++) {
		temp = c * x_data[ix] + s * x_data[iy];
		x_data[iy] = c * x_data[iy] - s * x_data[ix];
		x_data[ix] = temp;
		iy++;
		ix++;
	}
}

void xrotg(double *a, double *b, double *c, double *s)
{
	double roe;
	double absa;
	double absb;
	double scale;
	double ads;
	double bds;
	roe = *b;
	absa = fabs(*a);
	absb = fabs(*b);
	if (absa > absb) {
		roe = *a;
	}

	scale = absa + absb;
	if (scale == 0.0) {
		*s = 0.0;
		*c = 1.0;
		scale = 0.0;
		*b = 0.0;
	} else {
		ads = absa / scale;
		bds = absb / scale;
		scale *= sqrt(ads * ads + bds * bds);
		if (roe < 0.0) {
			scale = -scale;
		}

		*c = *a / scale;
		*s = *b / scale;
		if (absa > absb) {
			*b = *s;
		} else if (*c != 0.0) {
			*b = 1.0 / *c;
		} else {
			*b = 1.0;
		}
	}

	*a = scale;
}

void xscal(int n, double a, double x_data[], int ix0)
{
	int i9;
	int k;
	i9 = (ix0 + n) - 1;
	for (k = ix0; k <= i9; k++) {
		x_data[k - 1] *= a;
	}
}

void xswap(int n, double x_data[], int ix0, int iy0)
{
	int ix;
	int iy;
	int k;
	double temp;
	ix = ix0 - 1;
	iy = iy0 - 1;
	for (k = 1; k <= n; k++) {
		temp = x_data[ix];
		x_data[ix] = x_data[iy];
		x_data[iy] = temp;
		ix++;
		iy++;
	}
}
