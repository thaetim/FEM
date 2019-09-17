#pragma once

class Node
{
	double x, y, t;
	bool is_boundary;
public:
	Node();
	~Node();
	friend class Element;
	friend class Grid;
};

