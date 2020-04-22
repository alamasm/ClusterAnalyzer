#include "DBScanAlgorithm.h"
#include "algorithm"
void DBScanAlgorithm::find_clusters(vector<Point> points) {
    int n_clusters = 0;
    vector<int> label(points.size(), -1);//labels: 0 - noise, c - cluster

    for (size_t i = 0; i < points.size(); ++i) {
        if (label[i] != -1) {
            continue;
        }
        set<int> neighbors = get_neighbors(points, points[i], eps);
        if (neighbors.size() < min_pts) {
            label[i] = 0;
            continue;
        }
        n_clusters++;
        Cluster cluster;    
        label[i] = n_clusters;
        cluster.add_point(points[i]);
        for (size_t j : neighbors) {
            if (j == i) continue;
            
            if (label[j] == 0){
                 label[j] = n_clusters;
                 cluster.add_point(points[j]);
            }
            if (label[j] != -1){continue;}
            label[j] = n_clusters;
            cluster.add_point(points[j]);
            set<int> neighbors2 = get_neighbors(points, points[j], eps);
            if (neighbors2.size() >= min_pts){
                set_union(neighbors.begin(), neighbors.end(), neighbors2.begin(), neighbors2.end(), inserter(neighbors, neighbors.begin()));
            ;}
        }
        clusters.push_back(cluster);
    }
}

set<int> DBScanAlgorithm::get_neighbors(vector<Point> points, Point p, double eps) {
    set<int> neighbors;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].distance(p) < eps) neighbors.insert(i);
    }
    return neighbors;
}
