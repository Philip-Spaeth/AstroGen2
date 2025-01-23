#ifndef VEC4_H
#define VEC4_H

#include <iostream>

class vec4 {
public:
    float x, y, z, w;

    vec4();
    vec4(float x, float y, float z, float w);

    vec4(const vec4& v);
    vec4& operator=(const vec4& v);

    friend std::ostream& operator<<(std::ostream& os, const vec4& v);
};

#endif // VEC4_H
