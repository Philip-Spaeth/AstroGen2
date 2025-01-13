#pragma once

#include <iostream>
#include "vec3.h"

class mat4 {
public:
    // Konstruktoren
    mat4();
    mat4(float mat[4][4]);

    // Operatoren
    mat4 operator+(const mat4& other) const;
    mat4 operator-(const mat4& other) const;
    mat4 operator*(const mat4& other) const;
    float* operator*(const float vec[4]) const;

    // Methoden
    float determinant() const;
    mat4 inverse() const;
    mat4 transpose() const;

    static mat4 perspective(float fov, float aspect, float near, float far);
    static mat4 lookAt(const vec3& eye, const vec3& center, const vec3& up);

    // Hilfsfunktionen
    void print() const;
    const float* data() const;

private:
    // Matrix-Daten (4x4)
    float m[4][4];

    // Hilfsfunktionen zur Berechnung der Inversen
    void adjugate(mat4& adj) const;
    float minor(int row, int col) const;
    float cofactor(int row, int col) const;
};
