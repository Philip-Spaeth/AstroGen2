#include "mat4.h"
#include <cmath>
#include <iostream>

mat4::mat4() {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }
}
mat4::mat4(float mat[4][4]) {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m[i][j] = mat[i][j];
        }
    }
}

mat4 mat4::operator+(const mat4& other) const {
    mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = m[i][j] + other.m[i][j];
        }
    }
    return result;
}

mat4 mat4::operator-(const mat4& other) const {
    mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = m[i][j] - other.m[i][j];
        }
    }
    return result;
}

mat4 mat4::operator*(const mat4& other) const {
    mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = 0;
            for (int k = 0; k < 4; ++k) {
                result.m[i][j] += m[i][k] * other.m[k][j];
            }
        }
    }
    return result;
}

float* mat4::operator*(const float vec[4]) const {
    static float result[4];
    for (int i = 0; i < 4; ++i) {
        result[i] = 0;
        for (int j = 0; j < 4; ++j) {
            result[i] += m[i][j] * vec[j];
        }
    }
    return result;
}

float mat4::determinant() const {
    return m[0][0] * (
        m[1][1] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) -
        m[1][2] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) +
        m[1][3] * (m[2][1] * m[3][2] - m[2][2] * m[3][1])
    ) - m[0][1] * (
        m[1][0] * (m[2][2] * m[3][3] - m[2][3] * m[3][2]) -
        m[1][2] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) +
        m[1][3] * (m[2][0] * m[3][2] - m[2][2] * m[3][0])
    ) + m[0][2] * (
        m[1][0] * (m[2][1] * m[3][3] - m[2][3] * m[3][1]) -
        m[1][1] * (m[2][0] * m[3][3] - m[2][3] * m[3][0]) +
        m[1][3] * (m[2][0] * m[3][1] - m[2][1] * m[3][0])
    ) - m[0][3] * (
        m[1][0] * (m[2][1] * m[3][2] - m[2][2] * m[3][1]) -
        m[1][1] * (m[2][0] * m[3][2] - m[2][2] * m[3][0]) +
        m[1][2] * (m[2][0] * m[3][1] - m[2][1] * m[3][0])
    );
}

mat4 mat4::inverse() const {
    mat4 result;
    float det = determinant();
    if (det == 0) {
        std::cerr << "Matrix is singular and cannot be inverted" << std::endl;
        return result;
    }

    mat4 adj;
    adjugate(adj);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = adj.m[i][j] / det;
        }
    }

    return result;
}


void mat4::adjugate(mat4& adj) const {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            adj.m[j][i] = cofactor(i, j);
        }
    }
}

float mat4::minor(int row, int col) const {
    mat4 temp;
    int r = 0, c = 0;
    for (int i = 0; i < 4; ++i) {
        if (i == row) continue;
        c = 0;
        for (int j = 0; j < 4; ++j) {
            if (j == col) continue;
            temp.m[r][c++] = m[i][j];
        }
        ++r;
    }
    return temp.determinant();
}

float mat4::cofactor(int row, int col) const {
    return (row + col) % 2 == 0 ? minor(row, col) : -minor(row, col);
}

mat4 mat4::transpose() const {
    mat4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result.m[i][j] = m[j][i];
        }
    }
    return result;
}

mat4 mat4::perspective(float fov, float aspect, float near, float far) {
    mat4 result;
    float tanHalfFovy = tan(fov / 2.0f);

    result.m[0][0] = 1.0f / (aspect * tanHalfFovy);
    result.m[1][1] = 1.0f / tanHalfFovy;
    result.m[2][2] = -(far + near) / (far - near);
    result.m[2][3] = -1.0f;
    result.m[3][2] = -(2.0f * far * near) / (far - near);
    result.m[3][3] = 0.0f;

    return result;
}

mat4 mat4::lookAt(const vec3& eye, const vec3& center, const vec3& up) {
    vec3 f = (center - eye).normalize();
    vec3 s = f.cross(up).normalize();
    vec3 u = s.cross(f);

    mat4 result;
    result.m[0][0] = s.x;
    result.m[1][0] = s.y;
    result.m[2][0] = s.z;
    result.m[0][1] = u.x;
    result.m[1][1] = u.y;
    result.m[2][1] = u.z;
    result.m[0][2] = -f.x;
    result.m[1][2] = -f.y;
    result.m[2][2] = -f.z;
    result.m[3][0] = -s.dot(eye);
    result.m[3][1] = -u.dot(eye);
    result.m[3][2] = f.dot(eye);

    return result;
}

void mat4::print() const {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            std::cout << m[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

const float* mat4::data() const {
    return &m[0][0];
}
