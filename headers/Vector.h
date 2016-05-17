#ifndef VECTOR_H
#define VECTOR_H
#include <cmath>
#include <iostream>
class Vector
{
public:
    double x;
    double y;

    Vector();
    Vector(const Vector& obj);
    Vector(volatile const Vector& obj);
    Vector(double a, double b);
    Vector operator- () const;
    friend Vector operator+(const Vector& v1, const Vector& v2);
    friend Vector operator+(volatile const Vector& v1, volatile const Vector& v2);
    friend Vector operator-(const Vector& v1, const Vector& v2);
    friend Vector operator*(const Vector& v1, const Vector& v2);
    friend Vector operator*(double d, const Vector& v1);
    friend Vector operator*(const Vector& v1, double d);
    friend Vector operator/(const Vector& v1, double d);
    Vector& operator+=(const Vector &other);
    Vector& operator-=(const Vector &other);
    Vector& operator*=(const double d);
    double length() const;
    friend std::ostream& operator <<(std::ostream& os, Vector const& vector)
        {
            return os << "(" << vector.x << ", " << vector.y << ")";
        }
};

#endif /* VECTOR_H */

