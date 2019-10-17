#include "Group.h"
using namespace std;

Group::Group(int n, double x_center, double y_center, double x_dispersion, float y_dispersion) {
    this->x_center = x_center;
    this->y_center = y_center;
    this->x_dispersion = x_dispersion;
    this->y_dispersion = y_dispersion;
    vector<double> xs = get_random_normal_vector(n, -x_dispersion, x_dispersion);
    vector<double> ys = get_random_normal_vector(n, -y_dispersion, y_dispersion);
    for (int i = 0; i < n; ++i) {
        points.push_back(Point(x_center + xs[i], y_center + ys[n - i - 1]));
    }
}

void Group::rotate_relatively_to_center(double alpha) {
    double old_x;
    for (size_t i = 0; i < points.size(); ++i) {
        old_x = points[i].x;
        points[i].x = (old_x - x_center) * cos(alpha) - (points[i].y - y_center) * sin(alpha) + x_center;
        points[i].y = (old_x - x_center) * sin(alpha) + (points[i].y - y_center) * cos(alpha) + y_center;
    }
}

void Group::rotate_relatively_to_origin(double alpha) {
    double old_x;
    for (size_t i = 0; i < points.size(); ++i) {
        old_x = points[i].x;
        points[i].x = old_x * cos(alpha) - points[i].y * sin(alpha);
        points[i].y = old_x * sin(alpha) + points[i].y * cos(alpha);
    }
}


void Group::scale(double x_scale, double y_scale) {
    for (size_t i = 0; i < points.size(); ++i) {
        points[i].x = (points[i].x - x_center) * x_scale + x_center;
        points[i].y = (points[i].y - y_center) * y_scale + y_center;
    }
}

vector<double> Group::get_random_normal_vector(int n, double min, double max) {
    vector<double> res(n);
    uniform_real_distribution<double> unif(min,max);
    default_random_engine re;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 100; ++j) {
            res[i] += unif(re);
        }
        res[i] /= 100;
    }
    return res;
}