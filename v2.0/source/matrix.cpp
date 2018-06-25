#include "matrix.h"
#include <vector>
#include <cmath>
#include <string>

Matrix::~Matrix()
{
    this->n = 1;
    this->m = 1;
    this->a = {{0.0}};
}

Matrix::Matrix()
{
    this->n = 1;
    this->m = 1;
    this->a = {{0.0}};
}

Matrix::Matrix(int n1, int m1)
{
    this->n = n1;
    this->m = m1;
    this->a.resize(n1);
    for(int i = 0; i < n1; i++)
    {
        this->a[i].resize(m1, 0.0);
    }
}

Matrix::Matrix(std::vector<std::vector<double>> a1)
{
    this->n = a1.size();
    this->m = a1[0].size();
    this->a = a1;
}


void Matrix::operator =(Matrix B)
{
    this->n = B.n;
    this->m = B.m;
    this->a = B.a;
}

void Matrix::print()
{
    int N = this->n;
    int M = this->m;
    std::vector<std::vector<double>> A = this->a;
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < M; j++)
        {
            std::cout << A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    return;
}

Matrix matrixAdd(Matrix A, Matrix B)
{
    Matrix C = A;
    for(int i = 0; i < C.n; i++)
    {
        for(int j = 0; j < C.m; j++)
        {
            C.a[i][j] = C.a[i][j] + B.a[i][j];
        }
    }
    return C;
}

Matrix operator +(Matrix A, Matrix B)
{
    Matrix C = matrixAdd(A, B);
    return C;
}

Matrix operator -(Matrix A, Matrix B)
{
    Matrix C = (-1.0)*B;
    Matrix D = A + C;
    return D;
}

Matrix matrixProduct(Matrix A, Matrix B)
{
    Matrix C = Matrix(A.n, B.m);
    
    for(int i = 0; i < C.n; i++)
    {
        for(int j = 0; j < C.m; j++)
        {
            double x = 0.0;
            for(int k = 0; k < A.m; k++)
            {
                x = x + A.a[i][k]*B.a[k][j];
            }
            C.a[i][j] = x;
        }
    }
    
    return C;
}

Matrix operator *(Matrix A, Matrix B)
{
    Matrix C = matrixProduct(A, B);
    return C;
}

Matrix operator *(double p, Matrix A)
{
    Matrix B = scale(p, A);
    return B;
}

Matrix scale(double p, Matrix A)
{
    Matrix B = A;
    for(int i = 0; i < B.n; i++)
    {
        for(int j = 0; j < B.m; j++)
        {
            B.a[i][j] = p*B.a[i][j];
        }
    }
    return B;
}

Matrix zeroMatrix(int n, int m)
{
    Matrix A = Matrix(n, m);
    A = A - A;
    return A;
}

Matrix Id(int n)
{
    Matrix I = zeroMatrix(n, n);
    for(int i = 0; i < n; i++)
    {
        I.a[i][i] = 1.0;
    }
    return I;
}

Matrix RM(int n, int k, double p)
{
    k = k - 1;
    Matrix RM = Id(n);
    RM.a[k][k] = p;
    return RM;
}

Matrix RP(int n, int i, int j)
{
    i = i - 1;
    j = j - 1;
    Matrix RP = Id(n);
    std::vector<double> v = RP.a[i];
    RP.a[i] = RP.a[j];
    RP.a[j] = v;
    return RP;
}

Matrix RA(int n, int i, int j, double p)
{
    i = i - 1;
    j = j - 1;
    Matrix RA = Id(n);
    RA.a[i][j] = p;
    return RA;
}

Matrix Matrix::minor(int n1, int m1)
{
    n1 = n1 - 1;
    m1 = m1 - 1;
    std::vector<std::vector<double>> a1;
    a1.resize(this->n - 1);
    for(int i = 0; i < a1.size(); i++)
    {
        a1[i].resize(this->m - 1);
    }
    int i1 = 0;
    for(int i = 0; i < this->n; i++)
    {
        if(i != n1)
        {
            int j1 = 0;
            for(int j = 0; j < this->m; j++)
            {
                if(j != m1)
                {
                    a1[i1][j1] = this->a[i][j];
                    j1 = j1 + 1;
                }
            }
            i1 = i1 + 1;
        }
    }
    Matrix M = Matrix(a1);
    return M;
}

Matrix Matrix::transpose()
{
    Matrix B = *this;
    for(int i = 0; i < B.n; i++)
    {
        for(int j = 0; j < B.m; j++)
        {
            B.a[i][j] = this->a[j][i];
        }
    }
    return B;
}

double minorDet(Matrix A)
{
    double y = 0.0;
    if(A.m == 1) y = A.a[0][0];
    else
    {
        for(int k = 0; k < A.m; k++)
        {
            Matrix B = A.minor(1, k + 1);
            y = y + pow(-1.0, k)*A.a[0][k]*minorDet(B);
        }
    }
    return y;
}

Matrix cofactorMatrix(Matrix A)
{
    Matrix B = A;
    for(int i = 0; i < B.n; i++)
    {
        for(int j = 0; j < B.m; j++)
        {
            B.a[i][j] = pow(-1.0, i + j)*minorDet(A.minor(i + 1, j + 1));
        }
    }
    return B;
}

Matrix inverseMatrix(Matrix A, std::string method)
{
    Matrix B = A;
    if(method == "cf")
    {
        double d = minorDet(A);
        if(d == 0.0)
        {
            std::cout << "Can't inverse this matrix, sorry." << std::endl;
        }
        Matrix C = cofactorMatrix(A);
        Matrix adjoint = C.transpose();
        B = (1.0/d)*adjoint;
    } else if(method == "gj")
    {
        Matrix I = Id(A.n);
        for(int i = 0; i < A.n; i++)
        {
            Matrix T1;
            bool invertible = false;
            for(int k = i; k < A.n; k++)
            {
                if(B.a[k][i] != 0)
                {
                    T1 = RP(B.n, i + 1, k + 1);
                    k = B.n;
                    invertible = true;
                }
            }
            if(!invertible)
            {
                std::cout << "Can't invert this matrix, sorry." << std::endl;
                return B;
            }
            B = T1*B;
            I = T1*I;
            Matrix T2 = RM(B.n, i + 1, 1.0/B.a[i][i]);
            B = T2*B;
            I = T2*I;
            for(int k = 0; k < A.n; k++)
            {
                if(k != i){
                    Matrix Tk = RA(A.n, k + 1, i + 1, -B.a[k][i]);
                    B = Tk*B;
                    I = Tk*I;
                }
            }
        }
        B = I;
    } else
    {
        std::cout << "Method " << method << " not supported." << std::endl;
    }
    return B;
}
