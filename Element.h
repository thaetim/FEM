#pragma once
#include"Node.h"
#include"UniversalElement.h"
#include"GlobalData.h"
#include<iostream>

class Element
{
	int ID[4];
	double k, ro, c;
	double **H_el;
	double **C_el;
	double **H_BC_el;
	double *P_el;

  public:
	void calc_H_and_C_matrix(
		Node* array_of_nodes_for_this_element, 
		UniversalElement &universal_element, 
		int number_of_integration_points, 
		double* integration_points, 
		double* integration_weights);

	void calc_H_BC_and_P(
		Node* array_of_nodes_for_this_element,
		UniversalElement &universal_element);

	void add_H_BC_to_H();

	double** get_H_matrix();
	double** get_C_matrix();
	double** get_H_BC_matrix();
	double* get_P_vector();
	int* get_ID_arr();

	Element();
	~Element();
	friend class Grid;
};

