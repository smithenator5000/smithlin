#include "matrices.h"

using namespace std;

//---------------------------------------
Matrix::Matrix()
{
	destroy();
}
Matrix::Matrix(int n1, int m1)
{
	destroy();
	n = n1;
	m = m1;
	a = new double*[n];
	for(int i = 0; i < n; i++)
	{
		a[i] = new double[m];
		for(int j = 0; j < m; j++) a[i][j] = 0;
	}
}
Matrix::Matrix(int n1, int m1, double** a1, bool toDel)
{
	destroy();
	n = n1;
	m = m1;
	a = new double*[n];
	for(int i = 0; i < n; i++)
	{
		a[i] = new double[m];
		for(int j = 0; j < m; j++) a[i][j] = a1[i][j];
		if(toDel) delete [] a1[i];
	}
	if(toDel) 
	{
		for(int i = 0; i < n1; i++) delete [] a1[i];
		delete [] a1;
	}
}
void Matrix::operator =(Matrix B)
{
	destroy();
	n = B.n;
	m = B.m;
	a = new double*[n];
	for(int i = 0; i < n; i++)
	{
		a[i] = new double[m];
		for(int j = 0; j < m; j++) a[i][j] = B.a[i][j];
	}
	B.destroy();
	return;
}
void Matrix::operator /=(Matrix B)
{
	destroy();
	n = B.n;
	m = B.m;
	a = new double*[n];
	for(int i = 0; i < n; i++)
	{
		a[i] = new double[m];
		for(int j = 0; j < m; j++) a[i][j] = B.a[i][j];
	}
	return;
}
Matrix Matrix::operator +(Matrix B)
{
	Matrix C = Matrix(n, m);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++) C.a[i][j] = a[i][j] + B.a[i][j];
	}
	return C;
}
Matrix Matrix::operator -(Matrix B)
{
	Matrix C = Matrix(n, m);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++) C.a[i][j] = a[i][j] = B.a[i][j];
	}
	return C;
}
Matrix Matrix::operator *(Matrix B)
{
	Matrix C = Matrix(n, B.m);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < B.m; j++)
		{
			for(int k = 0; k < m; k++) C.a[i][j] += a[i][k]*B.a[k][j];
		}
	}
	return C;
}
Matrix Matrix::operator |(Matrix B)
{
	Matrix C = Matrix(n, m + B.m);
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++) C.a[i][j] = a[i][j];
		for(int j = 0; j < B.m; j++) C.a[i][m + j] = B.a[i][j];
	}
	return C;
}
Matrix Matrix::minor(int r, int s)
{
	Matrix C = Matrix(n - 1, m - 1);
	int ii = 0;
	for(int i = 0; i < n; i++)
	{
		if(i != r - 1)
		{
			for(int j = 0; j < s - 1; j++) C.a[ii][j] = a[i][j];
			for(int j = s; j < m; j++) C.a[ii][j - 1] = a[i][j];
			ii++;
		}
	}
	return C;
}
void Matrix::destroy()
{
	for(int i = 0; i < n; i++) delete [] a[i];
	if(n > 0) delete [] a;
	n = 0;
	m = 0;
}
void Matrix::print()
{
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < m; j++) cout << a[i][j] << "\t";
		cout << endl;
	}
}
//---------------------------------------
Matrix I(int n)
{
	Matrix I = Matrix(n, n);
	for(int i = 0; i < n; i++) I.a[i][i] = 1;
	return I;
}
Matrix RM(int n, int r, double c)
{
	Matrix RM = I(n);
	RM.a[r - 1][r - 1] = c;
	return RM;
}
Matrix RA(int n, int r, int s, double c)
{
	Matrix RA = I(n);
	RA.a[r - 1][s - 1] = c;
	return RA;
}
Matrix RP(int n, int r, int s)
{
	Matrix RP = I(n);
	double temp[n];
	for(int j = 0; j < n; j++)
	{
		temp[j] = RP.a[r - 1][j];
		RP.a[r - 1][j] = RP.a[s - 1][j];
		RP.a[s - 1][j] = temp[j];
	}
	return RP;
}
double i_det(Matrix A)
{
	double y;
	if(A.n == 1) y = A.a[0][0];
	else
	{
		y = 0;
		for(int j = 0; j < A.m; j++)
		{
			Matrix B = A.minor(1, j + 1);
			y += pow(-1.0, 2 + j)*A.a[0][j]*i_det(B);
			B.destroy();
		}
	}
	return y;
}
double gj_det(Matrix A1)
{
	double y;
	double k = 1;
	Matrix A;
	A /= A1;
	if(A.n == 1) y = A.a[0][0];
	else
	{
		int r;
		bool isTrivial = true;
		for(int i = 0; i < A.n; i++)
		{
			if(A.a[i][0] != 0) 
			{
				isTrivial = false;
				r = i;
			}
		}
		if(isTrivial) y = 0;
		else
		{
			k *= -1;
			Matrix P = RP(A.n, 1, r + 1);
			A = P*A;
			P.destroy();
			k *= 1.0/A.a[0][0];
			Matrix M = RM(A.n, 1, 1.0/A.a[0][0]);
			A = M*A;
			M.destroy();
			for(int i = 1; i < A.n; i++)
			{
				Matrix Ad = RA(A.n, i + 1, 1, -A.a[i][0]);
				A = Ad*A;
				Ad.destroy();
			}
			Matrix B = A.minor(1, 1);
			y = gj_det(B);
			B.destroy();
		}
	}
	y /= k;
	A.destroy();
	return y;
}
Matrix operator *(double k, Matrix A)
{
	Matrix C;
	C /= A;
	for(int i = 0; i < C.n; i++)
	{
		for(int j = 0; j < C.m; j++) C.a[i][j] *= k;
	}
	return C;
}
