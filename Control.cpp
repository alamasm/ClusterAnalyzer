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
    Algorithm* alg = new WaveAlgorithm(d);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_k_means(int k) {
    Algorithm* alg = new KMeansAlgorithm(k);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_spanning_tree(int k) {
    if (!spanning_tree_initialized) get_spanning_tree();
    Algorithm* alg = new SpanningTreeAlgorithm(k, &spanning_tree);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_hierarchical(int k) {
    Algorithm* alg = new HierarchicalAlgorithm(k);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_forel(double R) {
    Algorithm* alg = new FORELAlgorithm(R);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_dbscan(int min_pts, double eps) {
    Algorithm* alg = new DBScanAlgorithm(min_pts, eps);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
}

void Control::find_clusters_em(int k) {
    Algorithm* alg = new EMAlgorithm(k);
    alg->find_clusters(plane.get_points());
    finders.push_back(alg);
    em_index = finders.size() - 1;
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

pair<int, pair<vector<pair<Point, double>>, vector<Point>>> Control::get_em_data() {
    EMAlgorithm* em = (EMAlgorithm*) finders[em_index];
    return make_pair(em->k, make_pair(em->get_eighen(), em->mu));
}

vector<Cluster> Control::get_clusters(int n) {
    if (n == -1) return finders.back()->clusters;
    else return finders[n]->clusters;
}

int Control::set_groups(vector<Group> groups) {
    plane.groups = groups;
    return 1;
}
