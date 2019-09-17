#include"GlobalData.h"
#include"Element.h"
#include"Node.h"
#include"Grid.h"
#include"UniversalElement.h"
#include<iostream>
#include<cstdlib>

using namespace std;


void printSquareMatrix(double** M, int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			cout << M[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << endl;
}

void printVector(double* V, int size) {
	for (int i = 0; i < size; i++) {
		cout << V[i] << " ";
	}
	cout << endl << endl;
}

double* solveMatrixEquation(double** H, double* P, int dim) {

	double** A = new double*[dim];
	double* B = new double[dim];
	for (int i = 0; i < dim; i++) {
		A[i] = new double[dim];
		B[i] = P[i];
		for (int j = 0; j < dim; j++) {
			A[i][j] = H[j][i];
		}
	}

	for (int m = 0; m < dim; m++) {

		//find max value
		int max = m;
		for (int n = m + 1; n < dim; n++) {
			if (abs(A[n][m]) > abs(A[max][m])) {
				max = n;
			}
		}

		//swap row m with row with max in A and B
		/*double* tmp = A[m];
		A[m] = A[max];
		A[max] = tmp;*/
		swap(A[m], A[max]);
		swap(B[m], B[max]);

		if (abs(A[m][m]) <= 1E-10) {
			cout << "ERROR: MATRIX IS SINGULAR" << endl;
			return nullptr;
		}

		for (int i = m + 1; i < dim; i++) {
			double alpha = A[i][m] / A[m][m];
			B[i] -= alpha * B[m];
			for (int j = m; j < dim; j++) {
				A[i][j] -= alpha * A[m][j];
			}
		}
	}

	double* t1 = new double[dim];
	for (int i = dim - 1; i >= 0; i--) {
		double sum = 0.;
		for (int j = i + 1; j < dim; j++) {
			sum += A[i][j] * t1[j];
		}
		t1[i] = (B[i] - sum) / A[i][i];
	}

	cout << "equation solved" << endl;
	return t1;
}

double vectorMax(double* V, int dim) {
	int max = 0;
	for (int n = 1; n < dim; n++) {
		if (V[n] > V[max]) {
			max = n;
		}
	}
	return V[max];
}

double vectorMin(double* V, int dim) {
	int min = 0;
	for (int n = 1; n < dim; n++) {
		if (V[n] < V[min]) {
			min = n;
		}
	}
	return V[min];
}


int main() {

	//================ INITIALIZATION
	int a;
	GlobalData dat;
	Grid G(dat);
	UniversalElement UE;
	UE.calculateDiffs(dat.nP, dat.gauss_points);
	UE.calculateBoundaries(dat.H / (dat.nH-1), dat.L / (dat.nL-1), dat.nP, dat.gauss_points, dat.gauss_weights, dat.alfa, dat.temp_ambient);

	//check diffs calculated in UniversalElement
	UE.printDiffs(dat.nP);
	
	/*cout << endl << "Element 17: " << endl;
	G.Elements[16].calc_H_and_C_matrix(G.Nodes, UE, dat.nP, dat.gauss_points, dat.gauss_weights);
	cout << endl;

	cout << endl << "Element 1: " << endl;
	G.Elements[0].calc_H_BC_matrix(G.Nodes, UE, dat.nP, dat.gauss_points, dat.gauss_weights, dat.alfa);
	cout << endl;*/

	//================ STEADY STATE SOLUTION

	//number of elements
	int ne = (dat.nH - 1)*(dat.nL - 1);

	//calculate H, C and H_BC matrixes and P vector for each element
	for (int i = 0; i < ne; i++) {
		
		G.Elements[i].calc_H_and_C_matrix(G.Nodes, UE, dat.nP, dat.gauss_points, dat.gauss_weights);

		cout << "Element " << i << " H_BC:" << endl;
		G.Elements[i].calc_H_BC_and_P(G.Nodes, UE);

		if (i == 1) {
			std::cout << "Matrix H_el:" << endl;
			printSquareMatrix(G.Elements[i].get_H_matrix(), 4);
		}

		//add boundary conditions (H_BC) to H matrix -------------------------------------------- H_BC + H
		G.Elements[i].add_H_BC_to_H();
	}

	//print H C H_BC and P for the second element
	
	std::cout << "Matrix C_el:" << endl;
	printSquareMatrix(G.Elements[1].get_C_matrix(), 4);
	std::cout << "Matrix H_BC_el:" << endl;
	printSquareMatrix(G.Elements[1].get_H_BC_matrix(), 4);
	std::cout << "Vector P_el:" << endl;
	printVector(G.Elements[1].get_P_vector(), 4);
	

	//number of nodes
	int nn = dat.nH*dat.nL;

	//create global H and C matrixes and global P vector
	double** H = new double*[nn];
	double** C = new double*[nn];
	double* P = new double[nn];
	for (int i = 0; i < nn; i++) {
		H[i] = new double[nn];
		C[i] = new double[nn];
		P[i] = 0.;
		for (int j = 0; j < nn; j++) {
			H[i][j] = 0.;
			C[i][j] = 0.;
		}
	}

	//aggregate to global H and C matrixes
	double** H_el, **C_el, *P_el;
	int* ID_el;
	for (int i = 0; i < ne; i++) { // for each element
		ID_el = G.Elements[i].get_ID_arr();
		H_el = G.Elements[i].get_H_matrix();
		C_el = G.Elements[i].get_C_matrix();
		P_el = G.Elements[i].get_P_vector();
		for (int j = 0; j < 4; j++) {
			P[ID_el[j]] += P_el[j];
			for (int k = 0; k < 4; k++) {
				H[ID_el[j]][ID_el[k]] += H_el[j][k];
				C[ID_el[j]][ID_el[k]] += C_el[j][k];
			}
		}
	}

	//print global H matrix
	cout << endl << "GLOBAL H + H_BC MATRIX:" << endl;
	printSquareMatrix(H, nn);

	//print global C matrix
	cout << endl << "GLOBAL C MATRIX:" << endl;
	printSquareMatrix(C, nn);

	//print global P vector
	cout << endl << "GLOBAL P VECTOR:" << endl;
	printVector(P, nn);

	//================ TRANSIENT STATE SOLUTION

	int n_steps = dat.time / dat.time_step;
	double* t0_vec = new double[nn];
	double* t1_vec = new double[nn];
	for (int b = 0; b < nn; b++) {
		t0_vec[b] = dat.temp_initial;
	}

	//create [C]/dT and P_T = ([C]/dT)*{T0}
	double** C_T = new double*[nn];
	double* P_T = new double[nn];
	for (int i = 0; i < nn; i++) {
		C_T[i] = new double[nn];
		P_T[i] = 0.;
		for (int j = 0; j < nn; j++) {
			C_T[i][j] = 0.;
		}
	}

	int dT = dat.time_step;
	int T_sim = dat.time;

	for (int T = 0; T < T_sim; T += dT) {
		cout << endl << "TIME " << T << endl;

		// calculate [C]/dT
		for (int i = 0; i < nn; i++) {
			for (int j = 0; j < nn; j++) {
				C_T[i][j] = C[i][j] / dT;
			}
		}

		// [H]+[C]/dT and {P}+([C]/dT)*{T0}
		double t_tmp;
		for (int i = 0; i < nn; i++) {
			t_tmp = 0.;
			for (int j = 0; j < nn; j++) {
				if (T == 0) { H[i][j] += C_T[i][j]; }
				t_tmp += C_T[i][j] * t0_vec[j];
			}
			P_T[i] = P[i] + t_tmp;
		}

		//cout << T << endl;
		//print H matrix after 0 iteration
		if (T == 0) {
			cout << "[H]+[C]/dT:" << endl;
			printSquareMatrix(H,nn);

			cout << "{P}+([C]/dT)*{T0}:" << endl;
			printVector(P_T, nn);
		}

		//solve equation [H]{t1}+{P} -> find temperatures after dT
		t1_vec = solveMatrixEquation(H, P_T, nn);

		cout << "MIN AND MAX TEMP AFTER TIME " << T+dT << ": " 
			<< vectorMin(t1_vec, nn) << " " 
			<< vectorMax(t1_vec, nn) << endl << endl;

		//temperatures after this dT are the starting temperatures in the next dT
		for (int i = 0; i < nn; i++) {
			t0_vec[i] = t1_vec[i];
			t1_vec[i] = 0.;
		}
	}

	cout << endl;
	cin >> a;
	return 1;
}
