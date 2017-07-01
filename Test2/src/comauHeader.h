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

//real_T rtInf;
//real_T rtMinusInf;
//real_T rtNaN;
//real32_T rtInfF;
//real32_T rtMinusInfF;
//real32_T rtNaNF;

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
void AscissaCurvilineaOrientamento2(double T, const double phii[3], const double
		phif[3], double tf, emxArray_real_T *phie, emxArray_real_T *dphie,
		emxArray_real_T *ddphie);
void AscissaCurvilineaTraiettoria2(double T, const double qi[3], const double
		qf[3], double tf, emxArray_real_T *ps, emxArray_real_T *dps, emxArray_real_T
		*ddps);
void Click1PseudoinversaSamplesCpp_initialize();
void Click1PseudoinversaSamplesCpp_terminate();
double rt_roundd_snf(double u);
void Click1PseudoinversaSamplesCpp(const double q[6], const double dq[6], const
		double xd[6], const double dxd[6], double T, const char str[3], double K,
		double qf[6], double dqf_data[], int dqf_size[1]);
void invNxN(const double x[36], double y[36]);
void JacobianoGeometricoCpp(const double p[3], const double A[144], double Jg[36]);
double rt_atan2d_snf(double u0, double u1);
void kCpp(const double q[6], const char str[3], double p[3], double phi_data[],
		int phi_size[2], double R[9], double A[144]);
double norm(const double x[3]);
void eml_pinv(const double A_data[], const int A_size[2], double X_data[], int
		X_size[2]);
void power(const emxArray_real_T *a, emxArray_real_T *y);
void PVATraiettoria2(double T, double qf, double tf, emxArray_real_T *qd,
		emxArray_real_T *dqd, emxArray_real_T *ddqd);
void rt_InitInfAndNaN(size_t realSize);
boolean_T rtIsInf(real_T value);
boolean_T rtIsInfF(real32_T value);
boolean_T rtIsNaN(real_T value);
boolean_T rtIsNaNF(real32_T value);
real_T rtGetInf(void);
real32_T rtGetInfF(void);
real_T rtGetMinusInf(void);
real32_T rtGetMinusInfF(void);
real_T rtGetNaN(void);
real32_T rtGetNaNF(void);
void svd(const double A_data[], const int A_size[2], double U_data[], int
		U_size[2], double S_data[], int S_size[2], double V_data[], int V_size
		[2]);
emxArray_real_T *emxCreateND_real_T(int numDimensions, int *size);
emxArray_real_T *emxCreateWrapperND_real_T(double *data, int numDimensions, int *
		size);
emxArray_real_T *emxCreateWrapper_real_T(double *data, int rows, int cols);
emxArray_real_T *emxCreate_real_T(int rows, int cols);
void emxDestroyArray_real_T(emxArray_real_T *emxArray);
void emxInitArray_real_T(emxArray_real_T **pEmxArray, int numDimensions);
void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize);
void emxFree_real_T(emxArray_real_T **pEmxArray);
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
void Traj_initialize();
void Traj_terminate();
void Traj(double T, double tf, double tStop, const double p_i[3], const double
		p_f[3], const double phi_i[3], const double phi_f[3], emxArray_real_T *
		xd, emxArray_real_T *dxd, emxArray_real_T *ddxd);
void b_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0);
void c_xaxpy(int n, double a, const double x_data[], int ix0, double y_data[],
		int iy0);
void xaxpy(int n, double a, int ix0, double y_data[], int iy0);
double xdotc(int n, const double x_data[], int ix0, const double y_data[], int iy0);
double b_xnrm2(int n, const double x_data[], int ix0);
double xnrm2(int n, const double x_data[], int ix0);
void xrot(int n, double x_data[], int ix0, int iy0, double c, double s);
void xrotg(double *a, double *b, double *c, double *s);
void xscal(int n, double a, double x_data[], int ix0);
void xswap(int n, double x_data[], int ix0, int iy0);
