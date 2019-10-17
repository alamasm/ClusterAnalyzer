#include "Cluster.h"

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