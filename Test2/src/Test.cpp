//============================================================================
// Name        : Test.cpp
// Author      : Guarini - Iovine
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "comauHeader.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
using namespace std;
//
//void matTraj(emxArray_real_T *xd, emxArray_real_T *dxd, emxArray_real_T *ddxd, vector<vector<double> > matxd, vector<vector<double> > matdxd);
//void printfMatrix(vector<vector<double> > mat, int row, int col);

int main() {
	emxArray_real_T *xd;
	emxArray_real_T *dxd;
	emxArray_real_T *ddxd;

	double T=0.002;
	double tf= 5.0;
	double tStop=3.0;
	const double p_i[3]={-0.0000,4.0933,1.4835};
	const double p_f[3]={-0.0000, 3.3476, 1.0129};
	const double phi_i[3]={1.5708,-2.7053,0.0000};
	const double phi_f[3]={1.5708,-1.6581,0.0000};

	vector<vector<double> > matxd;
	vector<vector<double> > matdxd;

	emxInitArray_real_T(&xd, 2);
	emxInitArray_real_T(&dxd, 2);
	emxInitArray_real_T(&ddxd, 2);

	cout<<"\n"<<"GENERATION TRAJECTORY"<<endl;
	cout<<"LOADING ...\n"<<endl;
	Traj(T, tf, tStop, p_i, p_f, phi_i, phi_f, xd, dxd, ddxd);

	int row=xd->size[0];
	int col=xd->size[1];
	int index=0;
	vector<double> tmpxd;
	vector<double> tmpdxd;


	for(int c=0; c<col;c++){
		for(int r=0; r<row; r++){
			tmpxd.push_back(xd->data[index]);
			tmpdxd.push_back(dxd->data[index]);
			index++;
		}

		//tmpdxd=tmpxd;
		matxd.push_back(tmpxd);
		matdxd.push_back(tmpdxd);
		tmpxd.clear();
		tmpdxd.clear();
	}


	reverse(matxd.begin(),matxd.end());
	reverse(matdxd.begin(),matdxd.end());

	//	cout<<"\nPOSITION TRAJECTORY"<<endl;
	//	for(int i=0; i<col; i++){
	//		cout<<"COL: "<<col-i<<" -> ";
	//		for(int j=0; j<6; j++){
	//			cout<<matxd.at(i).at(j)<< " ";
	//		}
	//		cout<<endl;
	//
	//	}
	//
	//	cout<<"\nVELOCITY TRAJECTORY"<<endl;
	//	for(int i=0; i<col; i++){
	//		cout<<"COL: "<<i<<" -> ";
	//		for(int j=0; j<row; j++)
	//			cout<<" "<<matdxd[i][j]<<" ";
	//		cout<<endl;
	//	}

	/******************************************INVERSIONE***************************************************************************/
	//double matq[row][col];
	//double matdq[row][col];
	double q0[6]={1.5708, -2.4435, 0.6981, 0, 1.3090,  0};
	double dq0[6]={0, 0, 0, 0, 0, 0};
	int dqf_size[1];
	double dqf[6];
	double qf[6];
	const char str[3]={'Z','Y','Z'};
	double K=100;

	/****Vettori temporanei****/
	//double qtemp[row];
	//double dqtemp[row];
	double xdtemp[6];
	double dxdtemp[6];
	ofstream SaveFile("cpp_output.txt",std::ios_base::app);

	while(matxd.size()!=3601){
		//for(int i=0; i<row; i++){
		//	matq[i][0]=temp[i];
		//	matdq[i][0]=0;
		//}

		for(int i=0; i<6; i++){
			xdtemp[i]=matxd.at(matxd.size()-1).at(i);
			dxdtemp[i]=matdxd.at(matdxd.size()-1).at(i);
		}

		//cout<<"sizePrima: "<<matdxd.size()<<endl;

		//for(int k=1; k<col; k++){

		/*for(int i=0; i<row; i++){
			qtemp[i]=matq[i][k-1];
			dqtemp[i]=matdq[i][k-1];
			xdtemp[i]=matxd[i][k-1];
			dxdtemp[i]=matdxd[i][k-1];
		}*/

		Click1PseudoinversaSamplesCpp(q0, dq0, xdtemp, dxdtemp, T, str, K, qf, dqf, dqf_size);
		/*for(int i=0; i<row; i++){
			matq[i][k]=qf[i];
			matdq[i][k]=dqf[i];
		}*/


		cout<<endl;
		for(int i=0; i<6; i++){
			q0[i]=qf[i];
			dq0[i]=dqf[i];

			SaveFile<<q0[i]<<" ";

		}
		SaveFile<<"\n";

		//cout<<endl;
		matxd.pop_back();
		matdxd.pop_back();
		//cout<<"\nsizeDopo: "<<matdxd.size()<<endl;

	}
	//	cout<<"Inversione"<<endl;
	//	// STAMPA MATRICE TRAIETTORIA...
	//	for(int i=0; i<6; i++){
	//		cout<<"COLONNA: "<<i<<"";
	//		for(int j=0; j<row; j++)
	//			cout<<" "<<matq[j][i]<<" ";
	//		cout<<endl;
	//	}


	//}
	SaveFile.close();
	emxDestroyArray_real_T(ddxd);
	emxDestroyArray_real_T(dxd);
	emxDestroyArray_real_T(xd);

	return 0;
}

//void matTraj(emxArray_real_T *xd, emxArray_real_T *dxd, emxArray_real_T *ddxd, vector<vector<double> > matxd, vector<vector<double> > matdxd){
//	int index=0;
//	vector<double> tmpxd;
//	vector<double> tmpdxd;
//
//	for(int c=0; c<xd->size[1];c++){
//		for(int r=0; r<xd->size[0]; r++){
//			tmpxd.push_back(xd->data[index]);
//			tmpdxd.push_back(dxd->data[index]);
//			index++;
//		}
//		matxd.push_back(tmpxd);
//		matdxd.push_back(tmpdxd);
//	}
//	cout<<"DENTRO"<<endl;
//	cout<<matxd.size()<<endl;
//	cout<<matdxd.size()<<endl;
//}
//
//void printfMatrix(vector<vector<double> > mat, int row, int col){
//	for(int i=0; i<col; i++){
//		cout<<"COL: "<<i<<" -> ";
//		for(int j=0; j<row; j++)
//			cout<<" "<<endl;
//		//cout<<" "<<mat[j][i]<<" ";//vtmp.pop_back()<<" ";
//		cout<<endl;
//	}
//}
