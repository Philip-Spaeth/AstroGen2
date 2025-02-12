#include "vec4.h"

vec4::vec4() : x(0), y(0), z(0), w(0) {}

vec4::vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}

vec4::vec4(const vec4& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

vec4& vec4::operator=(const vec4& v) {
    if (this != &v) {
        x = v.x;
        y = v.y;
        z = v.z;
        w = v.w;
    }
    return *this;
}

std::ostream& operator<<(std::ostream& os, const vec4& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << ")";
    return os;
}