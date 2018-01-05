#ifndef VECTOR_H
#define VECTOR_H
#include "common.h"
class Vector
{
public:
	int n = 0; //dimension
	double* a; //elements
	Vector(); //default
	Vector(int); //zero vector
	Vector(int, double*, bool); //full definition
	Vector operator +(Vector);
	double operator *(Vector);
	Vector operator /(Vector);
	void operator =(Vector);
	void operator /=(Vector);
	void destroy(); //deallocates memory
	void print(); //print vector to console
};
Vector operator *(double, Vector);
#endif
