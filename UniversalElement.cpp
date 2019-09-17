#include "UniversalElement.h"
#include<iostream>


UniversalElement::UniversalElement()
{
	BC_down = new double*[2];
	BC_right = new double*[2];
	BC_up = new double*[2];
	BC_left = new double*[2];
	P_down = new double[2];
	P_right = new double[2];
	P_up = new double[2];
	P_left = new double[2];
	for (int i = 0; i < 2; i++) {
		BC_down[i] = new double[2];
		BC_right[i] = new double[2];
		BC_up[i] = new double[2];
		BC_left[i] = new double[2];
		P_down[i] = 0.;
		P_right[i] = 0.;
		P_up[i] = 0.;
		P_left[i] = 0.;
		for (int j = 0; j < 2; j++) {
			BC_down[i][j] = 0.;
			BC_right[i][j] = 0.;
			BC_up[i][j] = 0.;
			BC_left[i][j] = 0.;
		}
	}
}

void UniversalElement::calculateDiffs(unsigned int nP, double* points)
{
	unsigned int nP_squared = nP*nP;
	double ksi = 0., eta = 0.;
	N_values = new double*[nP_squared];
	N_ksi_diffs = new double*[nP_squared];
	N_eta_diffs = new double*[nP_squared];

	//arrays giving order of getting a point coordinate from gaussian quadrature
	int* ksi_order = new int[nP_squared];
	int* eta_order = new int[nP_squared];
	if (nP == 2) {
		ksi_order[0] = 0; ksi_order[1] = 1; ksi_order[2] = 1; ksi_order[3] = 0;
		eta_order[0] = 0; eta_order[1] = 0; eta_order[2] = 1; eta_order[3] = 1;
		/*int ksi_order[4] = { 0,1,1,0 };
		int eta_order[4] = { 0,0,1,1 };*/
	}
	else { //alteration for a 3 point variant
		ksi_order[0] = 0; ksi_order[1] = 1; ksi_order[2] = 2; ksi_order[3] = 2;
		ksi_order[4] = 2; ksi_order[5] = 1; ksi_order[6] = 0; ksi_order[7] = 0; ksi_order[8] = 1;
		eta_order[0] = 0; eta_order[1] = 0; eta_order[2] = 0; eta_order[3] = 1;
		eta_order[4] = 2; eta_order[5] = 2; eta_order[6] = 2; eta_order[7] = 1; eta_order[8] = 1;
		/*int ksi_order[9] = { 0,1,2,2,2,1,0,0,1 };
		int eta_order[9] = { 0,0,0,1,2,2,2,1,1 };*/
	}
	
	for (unsigned int i = 0; i < nP_squared; ++i) {

		ksi = points[ksi_order[i]];
		eta = points[eta_order[i]];

		N_values[i] = new double[4];
		N_values[i][0] = (1. - ksi)*(1. - eta) / 4.;
		N_values[i][1] = (1. + ksi)*(1. - eta) / 4.;
		N_values[i][2] = (1. + ksi)*(1. + eta) / 4.;
		N_values[i][3] = (1. - ksi)*(1. + eta) / 4.;

		N_ksi_diffs[i] = new double[4];
		N_ksi_diffs[i][0] = -(1. - eta) / 4.;
		N_ksi_diffs[i][1] = (1. - eta) / 4.;
		N_ksi_diffs[i][2] = (1. + eta) / 4.;
		N_ksi_diffs[i][3] = -(1. + eta) / 4.;

		N_eta_diffs[i] = new double[4];
		N_eta_diffs[i][0] = -(1. - ksi) / 4.;
		N_eta_diffs[i][1] = -(1. + ksi) / 4.;
		N_eta_diffs[i][2] = (1. + ksi) / 4.;
		N_eta_diffs[i][3] = (1. - ksi) / 4;
	}

	std::cout << "N_values, N_ksi_diffs and N_eta_diffs matrixes calculated" << endl << endl;
}

void UniversalElement::printDiffs(unsigned int nP)
{
	int n2 = nP*nP;

	cout << "N values:" << endl;
	for (int i = 0; i < n2; i++) {
		cout << "[	";
		for (int j = 0; j < 4; j++) {
			cout << N_values[i][j] << ",	";
		}
		cout << "]" << endl;
	}
	cout << endl;

	cout << "dN/dksi:" << endl;
	for (int i = 0; i < n2; i++) {
		cout << "[	";
		for (int j = 0; j < 4; j++) {
			cout << N_ksi_diffs[i][j] << ",	";
		}
		cout << "]" << endl;
	}
	cout << endl;

	cout << "dN/deta:" << endl;
	for (int i = 0; i < n2; i++) {
		cout << "[	";
		for (int j = 0; j < 4; j++) {
			cout << N_eta_diffs[i][j] << ",	";
		}
		cout << "]" << endl;
	}
	cout << endl;
}

void UniversalElement::calculateBoundaries(double el_height, double el_width, int nP, double* ipoints, double* iweights, double alfa, double temp_ambient) {
	
	double** BC_small = new double*[2];
	BC_small[0] = new double[2];
	BC_small[1] = new double[2];
	double* P_small = new double[2];

	for (int i = 0; i < 4; i++) {

		double l = 0.0;
		if (i==0 || i==2) {
			l = el_width;
		}
		else {
			l = el_height;
		}

		//empty the BC_small matrix for this boundary
		BC_small[0][0] = 0.;
		BC_small[0][1] = 0.;
		BC_small[1][0] = 0.;
		BC_small[1][1] = 0.;

		//empty the P_small vector for this boundary
		P_small[0] = 0.; P_small[1] = 0.;

		//integration on this boundary
		double det = l / 2; // dx/dksi = 1D Jakobian = ratio of natural (delta X) to universal (delta KSI) dimension
		for (int k = 0; k < nP; k++) { //for each integration point on this boundary
			double ksi = ipoints[k];
			double w = iweights[k];
			BC_small[0][0] += w*alfa*0.5*(1. - ksi)*0.5*(1. - ksi); //alfa*N1*N1
			double tmp = w*alfa*0.5*(1. - ksi)*0.5*(1. + ksi); //alfa*N2*N1=alfa*N2*N1
			BC_small[1][0] += tmp;
			BC_small[0][1] += tmp;
			BC_small[1][1] += w*alfa*0.5*(1. + ksi)*0.5*(1. + ksi); //alfa*N2*N2

			P_small[0] += w*alfa*temp_ambient*0.5*(1. - ksi);
			P_small[1] += w*alfa*temp_ambient*0.5*(1. + ksi);

			/*std::cout << "N1(ksi) = " << 0.5*(1 - ksi) << endl;
			std::cout << "N2(ksi) = " << 0.5*(1 + ksi) << endl;
			std::cout << "alfa = " << alfa << endl;*/
		}

		//multiply by ratio (divide by 1D Jakobian)
		BC_small[0][0] *= det;
		BC_small[0][1] *= det;
		BC_small[1][0] *= det;
		BC_small[1][1] *= det;
		P_small[0] *= det;
		P_small[1] *= det;

		//BC_small is now that small 2x2 matrix for this boundary
		//P_small is now that small size 2 vector for this boundary

		//--------------------BOUND MATRIXES PRINT
		std::cout << "bound " << i << ", det = " << det << std::endl << BC_small[0][0] << " " << BC_small[0][1] << std::endl << BC_small[1][0] << " " << BC_small[1][1] << std::endl << std::endl;
		
		//save the matrixes
		switch (i) {
		case 0:
			BC_down[0][0] = BC_small[0][0];
			BC_down[1][0] = BC_small[1][0];
			BC_down[0][1] = BC_small[0][1];
			BC_down[1][1] = BC_small[1][1];
			P_down[0] = P_small[0];
			P_down[1] = P_small[1];
			break;
		case 1:
			BC_right[0][0] = BC_small[0][0];
			BC_right[1][0] = BC_small[1][0];
			BC_right[0][1] = BC_small[0][1];
			BC_right[1][1] = BC_small[1][1];
			P_right[0] = P_small[0];
			P_right[1] = P_small[1];
			break;
		case 2:
			BC_up[0][0] = BC_small[0][0];
			BC_up[1][0] = BC_small[1][0];
			BC_up[0][1] = BC_small[0][1];
			BC_up[1][1] = BC_small[1][1];
			P_up[0] = P_small[0];
			P_up[1] = P_small[1];
			break;
		case 3:
			BC_left[0][0] = BC_small[0][0];
			BC_left[1][0] = BC_small[1][0];
			BC_left[0][1] = BC_small[0][1];
			BC_left[1][1] = BC_small[1][1];
			P_left[0] = P_small[0];
			P_left[1] = P_small[1];
			break;
		}
	}
	std::cout << "universal boundary matrixes and vectors calculated" << endl << endl;
}

double ** UniversalElement::getN_ksi_diffs()
{
	return N_ksi_diffs;
}

double ** UniversalElement::getN_eta_diffs()
{
	return N_eta_diffs;
}

double ** UniversalElement::getN_values()
{
	return N_values;
}

double ** UniversalElement::getH_BC(int id_of_boundary)
{
	switch (id_of_boundary) {
	case 0:
		return BC_down;
	case 1:
		return BC_right;
	case 2:
		return BC_up;
	case 3:
		return BC_left;
	}
	return nullptr;
}

double * UniversalElement::getP_vector(int id_of_boundary)
{
	switch (id_of_boundary) {
	case 0:
		return P_down;
	case 1:
		return P_right;
	case 2:
		return P_up;
	case 3:
		return P_left;
	}
	return nullptr;
}


UniversalElement::~UniversalElement()
{
}
