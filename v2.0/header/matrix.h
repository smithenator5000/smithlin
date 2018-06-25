#include <iostream>
#include <vector>
#include <string>

class Matrix
{
public:
    ~Matrix();
    Matrix();
    Matrix(int, int);
    Matrix(std::vector<std::vector<double>>);
    int n, m;
    std::vector<std::vector<double>> a;
    void print();
    void operator =(Matrix);
    Matrix minor(int, int);
    Matrix transpose();
};

Matrix matrixAdd(Matrix, Matrix);

Matrix operator +(Matrix, Matrix);

Matrix operator -(Matrix, Matrix);

Matrix matrixProduct(Matrix, Matrix);

Matrix operator *(Matrix, Matrix);

Matrix operator *(double, Matrix);

Matrix scale(double, Matrix);

Matrix zeroMatrix(int, int);

Matrix Id(int);

Matrix RM(int, int, double);

Matrix RP(int, int, int);

Matrix RA(int, int, int, double);

double minorDet(Matrix);

Matrix cofactorMatrix(Matrix);

Matrix inverseMatrix(Matrix, std::string method);