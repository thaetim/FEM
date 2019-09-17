#include "Element.h"



void Element::calc_H_and_C_matrix(Node* nodes, UniversalElement &ue, int nP, double* points, double* weights)
{
	double dxksi, dyksi, dxeta, dyeta, deter;
	double** Nkd = ue.getN_ksi_diffs();
	double** Ned = ue.getN_eta_diffs();
	double** N_vals = ue.getN_values();
	double **Nx_diffs, **Ny_diffs, **jakobians, **Nx_selfmatrix, **Ny_selfmatrix, **H_ip, **N_selfmatrix, **C_ip, *P_ip;
	int nP_sqr = nP*nP;
	Nx_diffs = new double*[nP_sqr];
	Ny_diffs = new double*[nP_sqr];
	jakobians = new double*[nP_sqr];
	Nx_selfmatrix = new double*[4];
	Ny_selfmatrix = new double*[4];
	N_selfmatrix = new double*[4];
	H_ip = new double*[4];
	C_ip = new double*[4];
	P_ip = new double[4];

	for (int ip = 0; ip < nP_sqr; ip++) {	//--------looping over integration points
		dxksi = 0.;
		dyksi = 0.;
		dxeta = 0.;
		dyeta = 0.;

		deter = 0.;

		//calc values in jakobi matrix and its determinant (DetJ)
		for (int l = 0; l < 4; l++) { 
			dxksi += Nkd[ip][l] * nodes[ID[l]].x;
			dyksi += Nkd[ip][l] * nodes[ID[l]].y;
			dxeta += Ned[ip][l] * nodes[ID[l]].x;
			dyeta += Ned[ip][l] * nodes[ID[l]].y;
		}
		deter = dxksi*dyeta - dyksi*dxeta;

		//-----------------JAKOBIAN DETERMINANT FOR THIS INTEGRATION POINT PRINT
		//std::cout << "Determinant for " << ip << " integration point: " << deter << endl;

		//calc the inverted jakobi matrix for this integration point (here as a row in the array)
		jakobians[ip] = new double[4];
		jakobians[ip][0] = dyeta / deter;
		jakobians[ip][1] = -dyksi / deter;
		jakobians[ip][2] = -dxeta / deter;
		jakobians[ip][3] = dxksi / deter;

		//----------------JAKOBI MATRIX FOR THIS INTEGRATION POINT PRINT
		/*
		std::cout << "Jakobian for " << ip + 1 << " integration point: " << 
			jakobians[ip][0] << "  " << jakobians[ip][1] << "	" << 
			jakobians[ip][2] << "  " << jakobians[ip][3] << std::endl;
		*/
		

		//calc the N/x and N/y differentials
		Nx_diffs[ip] = new double[4];
		Ny_diffs[ip] = new double[4];
		for (int j = 0; j < 4; j++) {
			Nx_diffs[ip][j] = Nkd[ip][j] * jakobians[ip][0] + Ned[ip][j] * (-jakobians[ip][1]);
			Ny_diffs[ip][j] = Nkd[ip][j] * (-jakobians[ip][2]) + Ned[ip][j] * jakobians[ip][3];
		}

		//calc the {dN/dx}{dN/dx}T and {dN/dy}{dN/dy}T -> self multiplied matrixes
		for (int i = 0; i < 4; i++) {
			Nx_selfmatrix[i] = new double[4];
			Ny_selfmatrix[i] = new double[4];
			N_selfmatrix[i] = new double[4];
			for (int j = 0; j < 4; j++) {
				Nx_selfmatrix[i][j] = Nx_diffs[ip][i] * Nx_diffs[ip][j];
				Ny_selfmatrix[i][j] = Ny_diffs[ip][i] * Ny_diffs[ip][j];
				N_selfmatrix[i][j] = N_vals[ip][i] * N_vals[ip][j];
			}
		}

		//arrays giving order of getting a point coordinate from gaussian quadrature
		int* ksi_order = new int[nP_sqr];
		int* eta_order = new int[nP_sqr];
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

		//calc the weight for this integration point
		double w = weights[ksi_order[ip]] * weights[eta_order[ip]];

		//---------INTEGRATION POINT COORDS AND WEIGHT PRINT
		/*
		std::cout << ip + 1 << " integration point: [" << points[ksi_order[ip]] << ", " << points[eta_order[ip]] <<
			"] with weight " << w << std::endl;
		*/

		//calc the w*k*( {dN/dx}{dN/dx}T + {dN/dy}{dN/dy}T )*DetJ -> local H matrix for this integration point ------------------------ WEIGHTS!!!! DODAæ
		//	AND add it to the general H matrix for this whole element
		for (int i = 0; i < 4; i++) {
			H_ip[i] = new double[4];
			C_ip[i] = new double[4];
			for (int j = 0; j < 4; j++) {
				H_ip[i][j] = w*k*(Nx_selfmatrix[i][j] + Ny_selfmatrix[i][j])*deter;
				H_el[i][j] += H_ip[i][j]; //adding to general H_el matrix

				C_ip[i][j] = w*c*ro*N_selfmatrix[i][j]*deter;
				C_el[i][j] += C_ip[i][j]; //adding to general C_el matrix
			}
		}

		//--------------------{dN/dx}{dN/dx}T AND {dN/dy}{dN/dy}T MATRIXES FOR THIS INTEGRATION POINT PRINT
		/*
		//print the self multiplied {dN/dx}{dN/dx}T matrix for this point
		std::cout << "{dN/dx}{dN/dx}T matrix for the " << ip+1 << " integration point:" << std::endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << Nx_selfmatrix[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		//print the self multiplied {dN/dy}{dN/dy}T matrix for this point
		std::cout << "{dN/dy}{dN/dy}T matrix for the " << ip+1 << " integration point:" << std::endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << Ny_selfmatrix[i][j] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
		*/


	}//--------end of looping over integration points

	//------------------------dN/dx DIFFERENTIALS FOR EACH INTEGRATION POINT PRINT
	/*std::cout << "dN/dx differentials for each integration point:" << std::endl;
	for (int i = 0; i < nP_sqr; i++) {
		std::cout << i + 1 << ": ";
		for (int j = 0; j < 4; j++) {
			std::cout << Nx_diffs[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout << "dN/dy differentials for each integration point:" << std::endl;
	for (int i = 0; i < nP_sqr; i++) {
		std::cout << i + 1 << ": ";
		for (int j = 0; j < 4; j++) {
			std::cout << Ny_diffs[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/

	//------------------------H MATRIX FOR THIS ELEMENT PRINT
	/*
	std::cout << "H matrix for this element:" << std::endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << H_el[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;*/

	//------------------------C MATRIX FOR THIS ELEMENT PRINT
	/*
	std::cout << "C matrix for this element:" << std::endl;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << C_el[i][j] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;*/

	//std::cout << "H_el and C_el matrixes calculated" << endl;
}

void Element::calc_H_BC_and_P(Node * nodes, UniversalElement & ue)
{
	//for every boundary (pair of nodes)
	for (int i = 0; i < 4; i++) {
		//here, indexes in A go like this: 0 1 2 3
		//here, indexes in B go like this: 1 2 3 0
		Node A = nodes[ID[i]];
		Node B = nodes[ID[(i + 1) % 4]];
		if (A.is_boundary && B.is_boundary) { //check if this pair is a boundary
			
			double** H_BC_tmp = ue.getH_BC(i);
			double* P_tmp = ue.getP_vector(i);
			switch (i) {
			case 0:
				//copy to upper left quarter
				for (int m = 0; m < 2; m++) {
					P_el[m] += P_tmp[m];
					for (int n = 0; n < 2; n++) {
						H_BC_el[m][n] += H_BC_tmp[m][n];
					}
				}
				std::cout << "down boundary copied" << endl;
				break;
			case 1:
				//copy to central four cells
				for (int m = 0; m < 2; m++) {
					P_el[m+1] += P_tmp[m];
					for (int n = 0; n < 2; n++) {
						H_BC_el[m+1][n+1] += H_BC_tmp[m][n];
					}
				}
				std::cout << "right boundary copied" << endl;
				break;
			case 2:
				//copy to lower right quarter
				for (int m = 0; m < 2; m++) {
					P_el[m+2] += P_tmp[m];
					for (int n = 0; n < 2; n++) {
						H_BC_el[m + 2][n + 2] += H_BC_tmp[m][n];
					}
				}
				std::cout << "up boundary copied" << endl;
				break;
			case 3:
				//copy to corner cells
				P_el[0] += P_tmp[0];
				P_el[3] += P_tmp[1];
				H_BC_el[0][0] += H_BC_tmp[0][0];
				H_BC_el[3][0] += H_BC_tmp[1][0];
				H_BC_el[0][3] += H_BC_tmp[0][1];
				H_BC_el[3][3] += H_BC_tmp[1][1];
				std::cout << "left boundary copied" << endl;
				break;
			}
		}
	}
	std::cout << endl;
	

	//std::cout << "H_BC_el matrix and P_el vector calculated" << endl;
}

void Element::add_H_BC_to_H() {

	//add H_BC_el to H_el
	for (int p = 0; p < 4; p++) {
		for (int r = 0; r < 4; r++) {
			H_el[p][r] += H_BC_el[p][r];
		}
	}
	//std::cout << "H_BC_el matrix and added to H_el matrix" << endl;
}

double ** Element::get_H_matrix()
{
	return H_el;
}

double ** Element::get_C_matrix()
{
	return C_el;
}

double ** Element::get_H_BC_matrix()
{
	return H_BC_el;
}

double * Element::get_P_vector()
{
	return P_el;
}

int * Element::get_ID_arr()
{
	return ID;
}

Element::Element()
{
	ID[0] = 0;
	ID[1] = 0;
	ID[2] = 0;
	ID[3] = 0;

	k = 1.;
	c = 1.;
	ro = 1.;

	H_el = new double*[4];
	C_el = new double*[4];
	H_BC_el = new double*[4];
	P_el = new double[4];
	for (int i = 0; i < 4; i++) {
		H_el[i] = new double[4];
		C_el[i] = new double[4];
		H_BC_el[i] = new double[4];
		P_el[i] = 0.;
		for (int j = 0; j < 4; j++) {
			H_el[i][j] = 0.;
			C_el[i][j] = 0.;
			H_BC_el[i][j] = 0.;
		}
	}
}

Element::~Element()
{
}
