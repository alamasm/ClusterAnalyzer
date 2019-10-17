#include "Plane.h"

using namespace std;

void Plane::add_group(int n, double x0, double y0, double x_dispersion, double y_dispersion) {
        Group g(n, x0, y0, x_dispersion, y_dispersion);
        groups.push_back(g);
}

int Plane::rotate_group_relatively_to_center(size_t group_n, double alpha) {
    if (group_n >= groups.size()) return -1;
    else {
        groups[group_n].rotate_relatively_to_center(alpha);
        return 1;
    }
}

int Plane::rotate_group_relatively_to_origin(size_t group_n, double alpha) {
    if (group_n >= groups.size()) return -1;
    else {
        groups[group_n].rotate_relatively_to_origin(alpha);
        return 1;
    }
}

vector<Point> Plane::get_points() {
    vector<Point> points;
    for (Group g : groups) {
        for (Point p : g.points) {
            points.push_back(p);
        }
    }
    return points;
}
