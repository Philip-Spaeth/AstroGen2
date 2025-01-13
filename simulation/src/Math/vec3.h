#pragma once

#include <iostream>
#include <cmath>

class vec3 {
public:
    double x, y, z;

    // Constructors
    vec3();
    vec3(double x, double y, double z);

    // Copy constructor
    vec3(const vec3& v);

    // Assignment operator
    vec3& operator=(const vec3& v);

    // Addition
    vec3 operator+(const vec3& v) const;
    vec3& operator+=(const vec3& v);

    // Subtraction
    vec3 operator-(const vec3& v) const;
    vec3& operator-=(const vec3& v);

    // Scalar multiplication
    vec3 operator*(double scalar) const;
    friend vec3 operator*(double scalar, const vec3& v);
    vec3& operator*=(double scalar);

    // Scalar division
    vec3 operator/(double scalar) const;
    vec3& operator/=(double scalar);

    // Length (magnitude) of the vector
    double length() const;

    //if statments
    bool operator==(const vec3& v) const;
    bool operator!=(const vec3& v) const;

    // Normalize the vector
    vec3 normalize() const;

    // Dot product
    double dot(const vec3& v) const;

    // Cross product
    vec3 cross(const vec3& v) const;

    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const vec3& v);

    // Convert to float array
    void toFloatArray(float arr[3]) const;

    const double* data() const { return &x; }
};
