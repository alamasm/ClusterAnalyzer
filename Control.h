#include <iostream>
#include <fstream>
#ifndef cluster_finder_include
#define cluster_finder_include
#include "ClusterFinder.h"
#endif

using namespace std;
class Control
{
    public:
    void create_group(int n, double x0, double y0, double x_dispersion, double y_dispersion);
    void rotate_group_relatively_to_center(int group_n, double alpha);
    void rotate_group_relatively_to_origin(int group_n, double alpha);
    void print_result(ofstream &out);
    void find_clusters(int algorithm, int d = 1);
    void print_clusters(ofstream &out);

    private:
    Plane plane;
    vector<Cluster> clusters;
    void print_points(vector<Point> points, ofstream &out, int group);
};