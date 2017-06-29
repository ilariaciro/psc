// PROTOTIPI METODI
/*void AscissaCurvilineaOrientamento2(double T, const double phii[3], const double
		phif[3], double tf, emxArray_real_T *phie, emxArray_real_T *dphie,
		emxArray_real_T *ddphie);
void AscissaCurvilineaTraiettoria2(double T, const double qi[3], const double
		qf[3], double tf, emxArray_real_T *ps, emxArray_real_T *dps, emxArray_real_T
		*ddps);*/
void Click1PseudoinversaSamplesCpp_initialize();
void Click1PseudoinversaSamplesCpp_terminate();
static double rt_roundd_snf(double u);
void Click1PseudoinversaSamplesCpp(const double q[6], const double dq[6], const
		double xd[6], const double dxd[6], double T, const char str[3], double K,
		double qf[6], double dqf_data[], int dqf_size[1]);
void invNxN(const double x[36], double y[36]);
void JacobianoGeometricoCpp(const double p[3], const double A[144], double Jg[36]);
static double rt_atan2d_snf(double u0, double u1);
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
