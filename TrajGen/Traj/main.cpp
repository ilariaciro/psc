//
// File: main.cpp
//
// MATLAB Coder version            : 3.0
// C/C++ source code generated on  : 28-Jun-2017 17:39:02
//

// Include Files
#include "rt_nonfinite.h"
#include "Traj.h"
#include "main.h"
#include "Traj_terminate.h"
#include "Traj_emxAPI.h"
#include "Traj_initialize.h"
#include <stdio.h>
#include <iostream>

// Function Declarations
static void argInit_3x1_real_T(double result[3]);
static double argInit_real_T();
static void main_Traj();
using namespace std;

// Function Definitions

//
// Arguments    : double result[3]
// Return Type  : void
//
static void argInit_3x1_real_T(double result[3])
{
	int idx0;

	// Loop over the array to initialize each element.
	for (idx0 = 0; idx0 < 3; idx0++) {
		// Set the value of the array element.
		// Change this value to the value that the application requires.
		result[idx0] = argInit_real_T();
	}
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
	return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_Traj()
{
	emxArray_real_T *xd;
	emxArray_real_T *dxd;
	emxArray_real_T *ddxd;
	double dv0[3];
	double dv1[3];
	double dv2[3];
	double dv3[3];
	emxInitArray_real_T(&xd, 2);
	emxInitArray_real_T(&dxd, 2);
	emxInitArray_real_T(&ddxd, 2);
	double T=0.002;
	double tf= 5.0;
	double tStop=3.0;
	const double p_i[3]={-0.0000,4.0933,1.4835};
	const double p_f[3]={-0.0000, 3.3476, 1.0129};
	const double phi_i[3]={1.5708,-2.7053,0.0000};
	const double phi_f[3]={1.5708,-1.6581,0.0000};

	cout<<"P1"<<endl;

	Traj(T, tf, tStop, p_i, p_f, phi_i, phi_f, xd, dxd, ddxd);

	cout<<"P3"<<endl;
	for(int i=0; i<xd->size[1]; i++){

	cout<<xd->data[i]<<endl;
	}

	emxDestroyArray_real_T(ddxd);
	emxDestroyArray_real_T(dxd);
	emxDestroyArray_real_T(xd);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
	// Initialize the application.
	// You do not need to do this more than one time.
	Traj_initialize();

	// Invoke the entry-point functions.
	// You can call entry-point functions multiple times.
	main_Traj();

	// Terminate the application.
	// You do not need to do this more than one time.
	Traj_terminate();
	return 0;
}

//
// File trailer for main.cpp
//
// [EOF]
//
