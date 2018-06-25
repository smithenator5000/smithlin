#ifndef MATRICES_H
#define MATRICES_H
#include "common.h"
class Matrix
{
public:
	int n = 0; //row number
	int m = 0; //column number
	double** a; //elements
	Matrix(); //default constructor
	Matrix(int, int); //zero matrix of definied dimensions
	Matrix(int, int, double**, bool); //full matrix definition
	void operator =(Matrix); //matrix equivalence
	void operator /=(Matrix); //matrix copy
	Matrix operator +(Matrix); //matrix addition
	Matrix operator -(Matrix); //matrix subtraction
	Matrix operator *(Matrix); //matrix multiplication
	Matrix operator |(Matrix); //matrix concatenation
	Matrix inverse(); //matrix inverse
	Matrix minor(int, int); //eliminates a defined row and column
	Matrix select(int, int, int, int); //isolates a section of a matrix defined by parameters
	void destroy(); //deallocates memory reserved for matrix
	void print(); //prints matrix to console
};
Matrix I(int); //n x n identity matrix
Matrix RM(int, int, double); //row multiplication matrix
Matrix RA(int, int, int, double); //row addition matrix
Matrix RP(int, int, int); //row permutation matrix
double i_det(Matrix); //compute determinant using iterative definition
double gj_det(Matrix); //compute determinant by simplifying problem with Gauss-Jordan elimination
Matrix operator*(double, Matrix);
#endif
