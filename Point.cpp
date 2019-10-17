#include <cmath>
#include "Point.h"

Point::Point(double x, double y) {
    this->x = x;
        this->y = y;
}

double Point::distance(Point p) {
    return sqrt((this->x - p.x) * (this->x - p.x) + (this->y - p.y) * (this->y - p.y));
}