#include <cmath>
#include "Point.h"

Point::Point(double x, double y) {
    this->x = x;
        this->y = y;
}

/*
Point::Point(const Point& point) {
    this->x = point.x;
    this->y = point.y;
}*/

Point::Point() {
    this->x = 0;
    this->y = 0;
}

double Point::distance(Point p) {
    return sqrt((this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y));
}

Point Point::operator+=(Point b) {
    this->x += b.x;
    this->y += b.y;
    return *this;
}

Point Point::operator/=(double x) {
    this->x /= x;
    this->y /= x;
    return *this;
}
