#include "Control.h"

using namespace std;

void Control::create_group(int n, double x0, double y0, double x_dispersion, double y_dispersion) {
    points_changed = 1;
    plane.add_group(n, x0, y0, x_dispersion, y_dispersion);
}

void Control::rotate_group_relatively_to_center(int group_n, double alpha) {
    points_changed = 1;
    plane.rotate_group_relatively_to_center(group_n, alpha);
}

void Control::rotate_group_relatively_to_origin(int group_n, double alpha) {
    points_changed = 1;
    plane.rotate_group_relatively_to_origin(group_n, alpha);
}

vector<Group> Control::get_groups() {
    return plane.groups;
}

/*
void Control::find_clusters(int algorithm, int d) {
    ClusterFinder finder = ClusterFinder(algorithm, d);
    clusters.push_back(finder.find_clusters(plane));
    clusters_searches_n++;
}*/

void Control::find_clusters_wave(double d) {
    ClusterFinder finder = ClusterFinder(0, d);
    finder.find_clusters_with_wave_algorithm(plane.get_points());
    finders.push_back(finder);
}

void Control::find_clusters_k_means(int k) {
    ClusterFinder finder = ClusterFinder(2);
    finder.find_clusters_with_k_means_algorithm(plane.get_points(), k);
    finders.push_back(finder);
}

void Control::find_clusters_spanning_tree(int k) {
    ClusterFinder finder = ClusterFinder(1);
    if (!spanning_tree_initialized) get_spanning_tree();
    finder.find_clusters_with_spanning_tree_algorithm(plane.get_points(), k, &spanning_tree);
    finders.push_back(finder);
}

void Control::find_clusters_hierarchical(int k) {
    ClusterFinder finder = ClusterFinder(3);
    finder.find_clusters_with_hierarchical_algorithm(plane.get_points(), k);
    finders.push_back(finder);
}

void Control::find_clusters_forel(double R) {
    ClusterFinder finder = ClusterFinder(4);
    finder.find_clusters_with_forel_algorithm(plane.get_points(), R);
    finders.push_back(finder);
}

void Control::find_clusters_dbscan(int min_pts, double eps) {
    ClusterFinder finder = ClusterFinder(5);
    finder.find_clusters_with_dbscan_algorithm(plane.get_points(), min_pts, eps);
    finders.push_back(finder);
}

pair<vector<vector<double>>, vector<Point>> Control::get_spanning_tree() {
    if (!spanning_tree_initialized || points_changed) {
        ClusterFinder finder = ClusterFinder(0, 0);
        spanning_tree = finder.spanning_tree(plane.get_points());
        spanning_tree_initialized = 1;
        points_changed = 0;
    }
    return spanning_tree;
}

vector<Cluster> Control::get_clusters(int n) {
    if (n == -1) return finders.back().clusters;
    else return finders[n].clusters;
}

int Control::set_groups(vector<Group> groups) {
    plane.groups = groups;
    return 1;
}
