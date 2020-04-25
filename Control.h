#include <iostream>
#include <fstream>
#ifndef cluster_finder_include
#define cluster_finder_include
#include "Algorithms.h"
#endif

using namespace std;
class Control
{
    public:
    int clusters_searches_n = 0;
    int points_changed = 0;
    int em_index;
    int spanning_tree_initialized = 0;
    pair<vector<vector<double>>, vector<Point>> spanning_tree;

    
    void create_group(int n, double x0, double y0, double x_dispersion, double y_dispersion);
    void rotate_group_relatively_to_center(int group_n, double alpha);
    void rotate_group_relatively_to_origin(int group_n, double alpha);
    int set_groups(vector<Group> groups);
    vector<Group> get_groups();
    
    //void find_clusters(int algorithm, int d = 1);
    void find_clusters_wave(double d);
    void find_clusters_k_means(int k);
    void find_clusters_spanning_tree(int k);
    void find_clusters_hierarchical(int k);
    void find_clusters_forel(double R);
    void find_clusters_dbscan(int min_pts, double eps);
    void find_clusters_em(int k);
    pair<vector<vector<double>>, vector<Point>> get_spanning_tree();
    pair<int, EM_step_data> get_em_data();
    pair<int, vector<EM_step_data>> get_em_animation_data();
    
    vector<Cluster> get_clusters(int n = -1);
    
    void print_clusters(ofstream &out);
    void print_result(ofstream &out);

    private:
    Plane plane;
    vector<Algorithm*> finders;
    void print_points(vector<Point> points, ofstream &out, int group);
};
