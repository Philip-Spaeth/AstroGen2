#ifndef VEC4_H
#define VEC4_H

#include <iostream>

class vec4 {
public:
    float x, y, z, w;

    // Constructors
    vec4();
    vec4(float x, float y, float z, float w);

    // Copy constructor
    vec4(const vec4& v);

    // Assignment operator
    vec4& operator=(const vec4& v);

    // Output stream operator
    friend std::ostream& operator<<(std::ostream& os, const vec4& v);
};

#endif // VEC4_H
