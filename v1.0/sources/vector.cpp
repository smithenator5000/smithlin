#include "vector.h"

using namespace std;

Vector::Vector()
{
	destroy();
}
Vector::Vector(int n1)
{
	destroy();
	n = n1;
	a = new double[n];
	for(int i = 0; i < n; i++) a[i] = 0;
}
Vector::Vector(int n1, double* a1, bool toDel)
{
	destroy();
	n = n1;
	a = new double[n];
	for(int i = 0; i < n; i++) a[i] = a1[i];
	if(toDel) delete [] a1;
}
Vector Vector::operator +(Vector v)
{
	Vector w = Vector(n);
	for(int i = 0; i < n; i++) w.a[i] = a[i] + v.a[i];
	return w;
}
Vector Vector::operator -(Vector v)
{
	Vector w = Vector(n);
	for(int i = 0; i < n; i++) w.a[i] = a[i] - v.a[i];
	return w;
}
void Vector::operator =(Vector v)
{
	destroy();
	n = v.n;
	a = new double[n];
	for(int i = 0; i < n; i++) a[i] = v.a[i];
	v.destroy();
	return;
}
void Vector::operator /=(Vector v)
{
	destroy();
	n = v.n;
	a = new double[n];
	for(int i = 0; i < n; i++) a[i] = v.a[i];
	return;
}
double Vector::operator *(Vector v)
{
	double y = 0;
	for(int i = 0; i < n; i++) y += a[i]*v.a[i];
	return y;
}
Vector Vector::operator <<(Vector v)
{
	Vector w = Vector(n);
	for(int i = 0; i < n; i++) w.a[i] = a[i] + v.a[i];
	v.destroy();
	destroy();
	return w;
}
Vector Vector::operator >>(Vector v)
{
	Vector w = Vector(n);
	for(int i = 0; i < n; i++) w.a[i] = a[i] - v.a[i];
	v.destroy();
	destroy();
	return w;
}
Vector Vector::operator /(Vector v)
{
	Vector w = Vector(3);
	w.a[0] = a[1]*v.a[2] - a[2]*v.a[1];
	w.a[1] = a[2]*v.a[0] - a[0]*v.a[2];
	w.a[2] = a[0]*v.a[1] - a[1]*v.a[0];
	return v;
}
double Vector::modulus()
{
	double y = 0;
	for(int i = 0; i < n; i++) y += a[i]*a[i];
	y = sqrt(y);
	return y;
}
void Vector::destroy()
{
	if(n > 0) delete [] a;
	n = 0;
}
Vector operator *(double k, Vector v)
{
	Vector w;
	w /= v;
	for(int i = 0; i < w.n; i++) w.a[i] *= k;
	return w;
}
