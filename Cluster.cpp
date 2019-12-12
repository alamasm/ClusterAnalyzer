#include "Cluster.h"
#include <math.h>
#include <iostream>
using namespace std;

Cluster::Cluster() {
    
}
Cluster::Cluster(vector<Point> points) {
    this->points = points;
} 

int Cluster::add_point(Point point) {
    this->points.push_back(point);
    return 1;
}

pair<Point, pair<Point, Point>> Cluster::get_eighen_vectors() {//center, v1, v2
    double a = 0;
    double b = 0;
    double c = 0;

    double xs = 0;
    double ys = 0;
    for (size_t i = 0; i < points.size(); ++i) {
        xs += points[i].x;
        ys += points[i].y;
    }
    double x_center = xs / points.size();
    double y_center = ys / points.size();
    for (size_t i = 0; i < points.size(); ++i) {
        a += (points[i].x - x_center) * (points[i].x - x_center);
        b += (points[i].x - x_center) * (points[i].y - y_center);
        c += (points[i].y - y_center) * (points[i].y - y_center);
    }
    a /= points.size();
    b /= points.size();
    c /= points.size();

    double l1 = 0.5 * (a + c + sqrt((a + c) * (a + c) - 4 * (a * c - b * b)));
    double l2 = 0.5 * (a + c - sqrt((a + c) * (a + c) - 4 * (a * c - b * b)));
    
    double x1 = b / (l1 - a);
    double y1 = 1 / (sqrt(x1 * x1 + 1)) * sqrt(2 * sqrt(l1));
    double x2 = b / (l2 - a);
    double y2 = 1 / (sqrt(x2 * x2 + 1)) * sqrt(2 * sqrt(l2));
    x1 *= y1;
    x2 *= y2;

    Point v1 = Point(x1, y1);
    Point v2 = Point(x2, y2);
    return make_pair(Point(x_center, y_center), make_pair(v1, v2));
}