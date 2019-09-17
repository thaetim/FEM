#include "GlobalData.h"
#include<iostream>
#include<cmath>

GlobalData::GlobalData()
{
	f.open("data.txt", ios::in);
	if (f.good() == true) {

		cout << "Reading data from file." << endl;
		string tmp;
		int index;

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		H = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		L = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		nH = (unsigned int)stoi(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		nL = (unsigned int)stoi(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		nP = (unsigned int)stoi(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		k = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		alfa = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		c = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		ro = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		temp_ambient = stod(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		time_step = stoi(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		time = stoi(tmp.substr(++index));

		getline(f, tmp);
		index = tmp.find_last_of(' ');
		temp_initial = stod(tmp.substr(++index));

		cout << endl
			<< "Global Data: " << endl
			<< "	H=" << H << endl 
			<< "	L=" << L << endl
			<< "	nH=" << nH << endl 
			<< "	nL=" << nL << endl
			<< "	nP=" << nP << endl
			<< "	k=" << k << endl
			<< "	alfa=" << alfa << endl
			<< "	c=" << c << endl
			<< "	ro=" << ro << endl
			<< "	temp_ambient=" << temp_ambient << endl
			<< "	time_step=" << time_step << endl
			<< "	time=" << time << endl
			<< "	temp_initial=" << temp_initial << endl;
		
		f.close();

		if(nP==2){
			gauss_points = new double[2];
			gauss_points[0] = -1. / sqrt(3.);
			gauss_points[1] = 1. / sqrt(3.);
			gauss_weights = new double[2];
			gauss_weights[0] = 1.;
			gauss_weights[1] = 1.;
			//std::cout << gauss_points[0] << ", " << gauss_weights[0];
		}else if(nP==3){
			gauss_points = new double[3];
			gauss_points[0] = -sqrt(3. / 5.);
			gauss_points[1] = 0.;
			gauss_points[2] = sqrt(3. / 5.);
			gauss_weights = new double[3];
			gauss_weights[0] = 5. / 9.;
			gauss_weights[1] = 0.;
			gauss_weights[2] = 5. / 9.;
			//std::cout << gauss_points[0] << ", " << gauss_weights[0];
		}else{
			cout << "GlobalData: ERROR: only 2 or 3 point Gauss quadrature methods are available"<<endl;
		}
	}
}


GlobalData::~GlobalData()
{
	delete [] gauss_points;
	delete [] gauss_weights;
}
