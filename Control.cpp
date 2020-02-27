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

void Control::find_clusters_wave(double d) {
    Algorithm alg = WaveAlgorithm(d);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_k_means(int k) {
    Algorithm alg = KMeansAlgorithm(k);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_spanning_tree(int k) {
    if (!spanning_tree_initialized) get_spanning_tree();
    Algorithm alg = SpanningTreeAlgorithm(k, &spanning_tree);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_hierarchical(int k) {
    Algorithm alg = HierarchicalAlgorithm(k);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_forel(double R) {
    Algorithm alg = FORELAlgorithm(R);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_dbscan(int min_pts, double eps) {
    Algorithm alg = DBScanAlgorithm(min_pts, eps);
    alg.find_clusters(plane.get_points());
    finders.push_back(alg);
}

pair<vector<vector<double>>, vector<Point>> Control::get_spanning_tree() {
    if (!spanning_tree_initialized || points_changed) {
        SpanningTreeAlgorithm alg = SpanningTreeAlgorithm(0, NULL);
        spanning_tree = alg.spanning_tree(plane.get_points());
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
