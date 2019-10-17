#include <vector>
#ifndef point_include
#define point_include
#include "Point.h"
#endif

class Cluster
{
    public:
    std::vector<Point> points;
    Cluster(std::vector<Point> points);
    Cluster();
    int add_point(Point point);
};