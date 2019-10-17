#include "Control.h"

using namespace std;

void Control::create_group(int n, double x0, double y0, double x_dispersion, double y_dispersion) {
    plane.add_group(n, x0, y0, x_dispersion, y_dispersion);
}

void Control::rotate_group_relatively_to_center(int group_n, double alpha) {
    plane.rotate_group_relatively_to_center(group_n, alpha);
}

void Control::rotate_group_relatively_to_origin(int group_n, double alpha) {
    plane.rotate_group_relatively_to_origin(group_n, alpha);
}

void Control::print_result(ofstream &out) {
    for (size_t i = 0; i < plane.groups.size(); ++i) {
        print_points(plane.groups[i].points, out, i);
    }
}

void Control::find_clusters(int algorithm, int d) {
    ClusterFinder finder = ClusterFinder(algorithm, d);
    clusters = finder.find_clusters(plane);
}

void Control::print_clusters(ofstream &out) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        print_points(clusters[i].points, out, i);
    }
}

void Control::print_points(vector<Point> points, ofstream &out, int group) {
    for (auto p : points) {
        out << p.x << " " << p.y << " " << group << endl;
    }
}