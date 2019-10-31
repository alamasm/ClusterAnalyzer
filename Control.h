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
    int clusters_searches_n = 0;
    int points_changed = 0;
    int spanning_tree_initialized = 0;
    pair<vector<vector<double>>, vector<Point>> spanning_tree;
    void create_group(int n, double x0, double y0, double x_dispersion, double y_dispersion);
    void rotate_group_relatively_to_center(int group_n, double alpha);
    void rotate_group_relatively_to_origin(int group_n, double alpha);
    void print_result(ofstream &out);
    void find_clusters(int algorithm, int d = 1);
    pair<vector<vector<double>>, vector<Point>> get_spanning_tree();
    void print_clusters(ofstream &out);
    vector<Cluster> get_clusters(int n = -1);
    vector<Group> get_groups();
    int set_groups(vector<Group> groups);

    private:
    Plane plane;
    vector<vector<Cluster>> clusters;
    void print_points(vector<Point> points, ofstream &out, int group);
};