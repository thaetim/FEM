#include "Grid.h"

Grid::Grid()
{
	
}

Grid::Grid(GlobalData & dat)
{
	nn = dat.nL*dat.nH;
	ne = (dat.nL-1)*(dat.nH-1);
	Elements = new Element[ne];
	Nodes = new Node[nn];

	//print numbers of nodes and elements
	std::cout << "Grid:" << endl
		<< "	nn=" << nn << endl
		<< "	ne=" << ne << endl << endl;

	//set coordinates for each node
	double dL = dat.L / (dat.nL - 1);
	double dH = dat.H / (dat.nH - 1);
	for (unsigned int i = 0; i < nn; i++) {
		Nodes[i].x = (i/dat.nH)*dL;
		Nodes[i].y = (i%dat.nH)*dH;
		cout << "Node " << i+1 << ": " << Nodes[i].x << ',' << Nodes[i].y << endl;
		//boundary condition
		if (Nodes[i].x == 0 || Nodes[i].x == dat.L || Nodes[i].y == 0 || Nodes[i].y == dat.H) {
			Nodes[i].is_boundary = true;
		}
	}

	//-------------BOUNDARY NODES PRINT
	/*for (unsigned int i = 0; i < nn; i++) {
		if (Nodes[i].x == 0 || Nodes[i].x == dat.L || Nodes[i].y == 0 || Nodes[i].y == dat.H) {
			std::cout << "Node " << i+1 << " is boundary" << endl;
		}
	}*/

	//cout << "Elements:" << endl;
	int first;
	for (unsigned int i = 0; i < ne; i++) {
		Elements[i].k = dat.k;
		Elements[i].c = dat.c;
		Elements[i].ro = dat.ro;

		first = i+i/(dat.nH-1);
		Elements[i].ID[0] = first;
		Elements[i].ID[1] = first + dat.nH;;
		Elements[i].ID[2] = first + dat.nH + 1;
		Elements[i].ID[3] = first + 1;

		//---------------------ELEMENT NODES COORDS PRINT
		/*
		if (i == 1) {
			cout << "Element " << i+1 << ':' << endl <<
				Elements[i].ID[0]+1 << ':' << '(' << Nodes[Elements[i].ID[0]].x << ',' << Nodes[Elements[i].ID[0]].y << ')' << endl <<
				Elements[i].ID[1]+1 << ':' << '(' << Nodes[Elements[i].ID[1]].x << ',' << Nodes[Elements[i].ID[1]].y << ')' << endl <<
				Elements[i].ID[2]+1 << ':' << '(' << Nodes[Elements[i].ID[2]].x << ',' << Nodes[Elements[i].ID[2]].y << ')' << endl << 
				Elements[i].ID[3]+1 << ':' << '(' << Nodes[Elements[i].ID[3]].x << ',' << Nodes[Elements[i].ID[3]].y << ')' << endl << endl;
		}
		*/
		//------------------------ELEMENT NODES IDS PRINT
		if (true) {
			cout << "Element " << i + 1 << ": [ " <<
				Elements[i].ID[0] + 1 << ' ' <<
				Elements[i].ID[1] + 1 << ' ' <<
				Elements[i].ID[2] + 1 << ' ' <<
				Elements[i].ID[3] + 1 << " ]" << endl;
		}
	}
}

Grid::~Grid()
{
}
