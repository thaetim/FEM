#pragma once
#include"GlobalData.h"
class UniversalElement
{
	double** N_ksi_diffs;
	double** N_eta_diffs;
	double** N_values;
	double** BC_down, **BC_right, **BC_up, **BC_left;
	double* P_down, *P_right, *P_up, *P_left;
public:
	UniversalElement();
	void calculateDiffs(unsigned int number_of_integration_points, double* gaussian_integration_points_tab);
	void printDiffs(unsigned int number_of_integration_points);
	void calculateBoundaries(
		double el_height, 
		double el_width, 
		int nP, 
		double* ipoints, 
		double* iweights, 
		double alfa,
		double temp_ambient);
	double** getN_ksi_diffs();
	double** getN_eta_diffs();
	double** getN_values();
	double** getH_BC(int id_orientation_of_boundary);
	double* getP_vector(int id_orientation_of_boundary);
	~UniversalElement();
	friend class Grid;
	friend class Element;
};

