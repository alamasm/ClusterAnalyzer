#include <vector>
using namespace std;
class Matrix {
    public:
    int n, m;
    Matrix(int n, int m, double c);
    bool set(int x, int y, double val);
    double get(int x, int y) const;
    vector<double>& operator[](size_t i) {return data[i];}
    Matrix dot(const Matrix& b) const;
    double det() const;
    Matrix operator*(double a);
    private:
    vector<vector<double>> data;
    bool check_indexes(int x, int y) const;
};