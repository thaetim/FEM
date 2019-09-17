#pragma once
#include<fstream>
#include<iostream>
#include<string>

using namespace std;

class GlobalData
{
public:
	double *gauss_points, *gauss_weights;
	fstream f;
	double H, L, k, alfa, c, ro, temp_ambient, temp_initial;
	unsigned int nH, nL, nP;
	int time_step, time;
	friend class Grid;
	friend class Element;
	friend class UniversalElement;
	GlobalData();
	~GlobalData();
};

