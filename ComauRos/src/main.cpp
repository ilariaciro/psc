// Include Files
/*#include "rt_nonfinite.h"
#include "Traj.h"
#include "main.h"
#include "Traj_terminate.h"
#include "Traj_emxAPI.h"
#include "Traj_initialize.h"
#include "Click1PseudoinversaSamplesCpp.h"
#include "main.h"
#include "Click1PseudoinversaSamplesCpp_terminate.h"
#include "Click1PseudoinversaSamplesCpp_initialize.h"




// Function Declarations
//static void argInit_3x1_real_T(double result[3]);
//static double argInit_real_T();

//static void main_Traj();

// Function Definitions

//
// Arguments    : double result[3]
// Return Type  : void
//

//static void argInit_3x1_real_T(double result[3])
//{
//  int idx0;
//
//  // Loop over the array to initialize each element.
//  for (idx0 = 0; idx0 < 3; idx0++) {
//    // Set the value of the array element.
//    // Change this value to the value that the application requires.
//    result[idx0] = argInit_real_T();
//  }
//}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}
 */



//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//

//void printMatrixTraj(double matTraj[][], int row, int col);

#include <stdio.h>
#include <iostream>

#include "comauHeader.h"
using namespace std;

int main(int, const char * const [])
{
	// Initialize the application.
	// Traj_initialize();
	//Click1PseudoinversaSamplesCpp_initialize();

	// Invoke the entry-point functions.
	// You can call entry-point functions multiple times.
	//main_Traj();
	emxArray_real_T *xd;
	emxArray_real_T *dxd;
	emxArray_real_T *ddxd;

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

	cout<<"\n"<<"GENERAZIONE TRAIETTORIA"<<endl;
	Traj(T, tf, tStop, p_i, p_f, phi_i, phi_f, xd, dxd, ddxd);

	// COSTRUZIONE MATRICE TRAIETTORIA
	int row=xd->size[0];
	int col=xd->size[1];
	double matxd[row][col];
	double matdxd[row][col];

	int index=0;
	for(int c=0; c<col; c++)
		for(int r=0; r<row; r++){
			matxd[r][c]=xd->data[index];
			matdxd[r][c]=dxd->data[index];
			index++;
		}

	/*	cout<<"Traiettoria"<<endl;
	// STAMPA MATRICE TRAIETTORIA...
	for(int i=0; i<col; i++){
		cout<<"COLONNA: "<<i<<"";
		for(int j=0; j<row; j++)
			cout<<" "<<matxd[j][i]<<" ";
		cout<<endl;
	}

	cout<<"Velocità"<<endl;
	// STAMPA MATRICE TRAIETTORIA...
	for(int i=0; i<col; i++){
		cout<<"COLONNA: "<<i<<"";
		for(int j=0; j<row; j++)
			cout<<" "<<matdxd[j][i]<<" ";
		cout<<endl;
	}*/

	/********************************* INVERSIONE *********************************************************************/
	double matq[row][col];
	double matdq[row][col];
	double temp[6]={1.5708, -2.4435, 0.6981, 0, 1.3090,  0};
	int dqf_size[1];
	double dqf[row];
	double qf[row];
	const char str[3]={'Z','Y','Z'};
	double K=100;

	/****Vettori temporanei****/
	double qtemp[row];
	double dqtemp[row];
	double xdtemp[row];
	double dxdtemp[row];

	for(int i=0; i<row; i++){
		matq[i][0]=temp[i];
		matdq[i][0]=0;
	}

	for(int k=1; k<col; k++){

		for(int i=0; i<row; i++){
			qtemp[i]=matq[i][k-1];
			dqtemp[i]=matdq[i][k-1];
			xdtemp[i]=matxd[i][k-1];
			dxdtemp[i]=matdxd[i][k-1];
		}
		Click1PseudoinversaSamplesCpp(qtemp, dqtemp, xdtemp, dxdtemp, T, str, K, qf, dqf, dqf_size);
		for(int i=0; i<row; i++){
			matq[i][k]=qf[i];
			matdq[i][k]=dqf[i];
		}

	}
	cout<<"Inversione"<<endl;
		// STAMPA MATRICE TRAIETTORIA...
		for(int i=0; i<col; i++){
			cout<<"COLONNA: "<<i<<"";
			for(int j=0; j<row; j++)
				cout<<" "<<matq[j][i]<<" ";
			cout<<endl;
		}



	emxDestroyArray_real_T(ddxd);
	emxDestroyArray_real_T(dxd);
	emxDestroyArray_real_T(xd);

	// Terminate the application.
	//Traj_terminate();
	//Click1PseudoinversaSamplesCpp_terminate();

	return 0;
}

/*void printMatrixTraj(double matTraj[][], int row, int col){
  for(int i=0; i<col; i++){
      for(int j=0; j<row; j++)
        cout<<" "<<matTraj[j][i]<<" ";
      cout<<endl;
    }
}
 */
//
// File trailer for main.cpp
//
// [EOF]
//
