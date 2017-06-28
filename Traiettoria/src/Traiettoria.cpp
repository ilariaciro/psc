//============================================================================
// Name        : Traiettoria.cpp
// Author      : Guarini - Iovine
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "Traj.h"
#include "Traj_emxAPI.h"
#include "Traj_emxutil.h"
#include "Traj_initialize.h"
#include "Traj_terminate.h"
#include "Traj_types.h"
#include <iostream>
using namespace std;

int main() {

	double T=0.002;
	double tf= 5.0;
	double tStop=3.0;
	const double p_i[3]={-0.0000,4.0933,1.4835};
	const double p_f[3]={-0.0000, 3.3476, 1.0129};
	const double phi_i[3]={1.5708,-2.7053,0.0000};
	const double phi_f[3]={1.5708,-1.6581,0.0000};
	emxArray_real_T *xd=NULL;
	emxArray_real_T *dxd=NULL;
	emxArray_real_T *ddxd=NULL;

	cout<<"\n"<<"Pen1"<<endl;
	Traj(T, tf, tStop, p_i, p_f, phi_i, phi_f, xd, dxd, ddxd);
	cout<<"\n"<<"Pen"<<endl;
	cout<<"\n"<<"QUI"<<endl;
	return 0;
}
