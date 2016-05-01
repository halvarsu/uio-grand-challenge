#include "headers/Vector.h"

Vector::Vector(): x(0.0), y(0.0)
{}

Vector::Vector(const Vector& obj): x(obj.x), y(obj.y)
{}

Vector::Vector(double a, double b): x(a), y(b)
{}

Vector Vector::operator- () const
{
    return Vector(-x, -y);
}

Vector operator+(const Vector& v1, const Vector& v2)
{
    return Vector(v1.x+v2.x, v1.y+v2.y);
}

Vector operator-(const Vector& v1, const Vector& v2)
{
    return Vector(v1.x - v2.x, v1.y - v2.y);
}

Vector operator*(const Vector& v1, const Vector& v2)
{
    return Vector(v1.x*v2.x, v1.y*v2.y);
}

Vector operator*(double d, const Vector& v1)
{
    return Vector(v1.x*d, v1.y*d);
}

Vector operator*(const Vector& v1, double d){
    return Vector(v1.x*d, v1.y*d);
}

Vector operator/(const Vector& v1, double d)
{
    if(d != 0)
        return Vector(v1.x/d, v1.x/d);
    else
        return Vector(NAN,NAN); // Poor solution
}

Vector& Vector::operator+=(const Vector &other)
{
    this->x += other.x;
    this->y += other.y;
    return *this;
}

Vector& Vector::operator-=(const Vector &other)
{
    this->x -= other.x;
    this->y -= other.y;
    return *this;
}

double Vector::length() const
{
    if(x==0.0)
        return y;
    else if(y==0.0)
        return x;
    else
        return sqrt(x*x + y*y);
}
