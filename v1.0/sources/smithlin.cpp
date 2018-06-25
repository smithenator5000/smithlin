#include "common.h"
#include "smithlin.h"

using namespace std;

Vector operator *(Matrix A, Vector v)
{
	Vector w = Vector(A.n);
	for(int i = 0; i < A.n; i++)
	{
			for(int k = 0; k < A.m; k++) v.a[i] += A.a[i][k]*v.a[k];
	}
	return w;
}
Matrix operator |(Matrix A, Vector v)
{
	Matrix B = Matrix(A.n, A.m + 1);
	for(int i = 0; i < B.n; i++)
	{
		for(int j = 0; j < A.m; j++) B.a[i][j] = A.a[i][j];
		B.a[i][A.m] = v.a[i];
	}
	return B;
}
Matrix operator |(Matrix A, Matrix B)
{
	Matrix C = Matrix(A.n, A.m + B.m);
	for(int i = 0; i < C.n; i++)
	{
		for(int j = 0; j < A.m; j++) C.a[i][j] = A.a[i][j];
		for(int j = 0; j < B.m; j++) C.a[i][j + A.m] = B.a[i][j];
	}
	return C;
}
bool isBasis(int n, Vector* B, bool toDel)
{
	Matrix A = Matrix(n, 0);
	for(int i = 0; i < n; i++) A.a[i][0] = B[0].a[i];
	for(int i = 1; i < n; i++)
	{
		A = A|B[i];
	}
	double y = gj_det(A);
	A.destroy();
	if(y == 0) return false;
	else return true;
}
using namespace std;
