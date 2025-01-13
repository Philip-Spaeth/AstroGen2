#include "vec3.h"
#include <iostream>
#include <cmath>

// Default constructor
vec3::vec3() : x(0), y(0), z(0) {}

// Parameterized constructor
vec3::vec3(double x, double y, double z) : x(x), y(y), z(z) {}

// Copy constructor
vec3::vec3(const vec3& v) : x(v.x), y(v.y), z(v.z) {}

// Assignment operator
vec3& vec3::operator=(const vec3& v) {
    if (this == &v) return *this;
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
}

// Addition
vec3 vec3::operator+(const vec3& v) const {
    return vec3(x + v.x, y + v.y, z + v.z);
}

vec3& vec3::operator+=(const vec3& v) {
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

// Subtraction
vec3 vec3::operator-(const vec3& v) const {
    return vec3(x - v.x, y - v.y, z - v.z);
}

vec3& vec3::operator-=(const vec3& v) {
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

// Scalar multiplication
vec3 vec3::operator*(double scalar) const {
    return vec3(x * scalar, y * scalar, z * scalar);
}

vec3 operator*(double scalar, const vec3& v) {
    return vec3(v.x * scalar, v.y * scalar, v.z * scalar);
}

vec3& vec3::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
}

// Scalar division
vec3 vec3::operator/(double scalar) const {
    return vec3(x / scalar, y / scalar, z / scalar);
}

vec3& vec3::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

// Length (magnitude) of the vector
double vec3::length() const {
    return std::sqrt(x * x + y * y + z * z);
}

// Normalize the vector
vec3 vec3::normalize() const 
{
    double len = length();
    if (len > 0) {
        return *this / len;
    }
    return vec3(0, 0, 0);
}

// Dot product
double vec3::dot(const vec3& v) const {
    return x * v.x + y * v.y + z * v.z;
}

// Cross product
vec3 vec3::cross(const vec3& v) const {
	return vec3(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
}

//if statments
bool vec3::operator==(const vec3& v) const {
    return x == v.x && y == v.y && z == v.z;
}

bool vec3::operator!=(const vec3& v) const {
    return x != v.x || y != v.y || z != v.z;
}

// Output stream operator
std::ostream& operator<<(std::ostream& os, const vec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

// Convert to float array
void vec3::toFloatArray(float arr[3]) const {
    arr[0] = static_cast<float>(x);
    arr[1] = static_cast<float>(y);
    arr[2] = static_cast<float>(z);
}
