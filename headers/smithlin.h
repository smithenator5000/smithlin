#ifndef SMITHLIN_H
#define SMITHLIN_H
#include "matrices.h"
#include "vector.h"
Vector operator *(Matrix, Vector);
Matrix operator |(Matrix, Vector);
Matrix operator |(Vector, Vector);
bool isBasis(int, Vector*, bool);
#endif
