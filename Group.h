#include <vector>
#include <cmath>
#include <random>
#ifndef point_include
#define point_include
#include "Point.h"
#endif
using namespace std;
class Group
{
    public:
    vector<Point> points;
    double x_center;
    double y_center;
    double x_dispersion;
    double y_dispersion;
    Group(int n, double x_center, double y_center, double x_dispersion, float y_dispersion);
    void rotate_relatively_to_center(double alpha);
    void rotate_relatively_to_origin(double alpha);
    void scale(double x_scale, double y_scale);
    private:
    vector<double> get_random_normal_vector(int n, double min, double max);
};