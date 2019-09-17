#pragma once
#include"GlobalData.h"
#include"Element.h"
#include"Node.h"
#include"UniversalElement.h"

#include<iostream>

class Grid
{
public:
	unsigned int ne, nn;
	Element* Elements;
	Node* Nodes;
	Grid();
	Grid(GlobalData &dat);
	~Grid();
	friend class Element;
};

