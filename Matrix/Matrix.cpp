#include "Matrix.h"

Matrix::Matrix(int n, int m, double c) {
    this->n = n;
    this->m = m;
    for (int i = 0; i < n; ++i) {
        data.push_back(vector<double>(m));
        for (int j = 0; j < m; ++j) {
            data[i][j] = c;
        }
    }
}

Matrix Matrix::operator*(double a) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            data[i][j] *= a;
        }
    }
    return *this;
}

double Matrix::get(int x, int y) const {
    if (check_indexes(x, y)) return data[x][y];
    else return -1e9;
}

bool Matrix::set(int x, int y, double val) {
    if (!check_indexes(x, y)) return 0;
    data[x][y] = val;
    return 1;
}

bool Matrix::check_indexes(int x, int y) const {
    return (x >= 0 && x < n && y >= 0 && y < m);
}

Matrix Matrix::dot(const Matrix& b) const {
    if (m != b.n) return Matrix(0, 0, 0);
    Matrix res = Matrix(n, b.m, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < b.m; ++j) {
            double c = 0;
            for (int k = 0; k < m; ++k) {
                c += get(i, k) * b.get(k, j);
            }
            res[i][j] = c;
        }
    }
    return res;
}