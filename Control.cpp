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

void Control::find_clusters(int algorithm, int d) {
    ClusterFinder finder = ClusterFinder(algorithm, d);
    clusters.push_back(finder.find_clusters(plane));
    clusters_searches_n++;
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
    if (n == -1) return clusters.back();
    else return clusters[n];
}

int Control::set_groups(vector<Group> groups) {
    plane.groups = groups;
    return 1;
}