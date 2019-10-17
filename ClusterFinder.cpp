#include "ClusterFinder.h"

ClusterFinder::ClusterFinder(int algorithm, int d) {
    this->algorithm = algorithm;
    this->d = d;
}

vector<Cluster> ClusterFinder::find_clusters(Plane plane) {
    switch (algorithm) {
        case ALGORITHM_WAVE:
            return find_clusters_with_wave_algorithm(plane.get_points());
        case ALGORITHM_SPANNING_TREE:
            return find_clusters_with_spanning_tree_algorithm(plane.get_points());
        case ALGORITHM_K_MEANS:
            return find_clusters_with_k_means_algorithm(plane.get_points());
    }
    return vector<Cluster>();
}

vector<Cluster> ClusterFinder::find_clusters_with_wave_algorithm(vector<Point> points) {
    vector<vector<int>> graph(points.size());
    vector<Cluster> clusters;
    for (size_t i = 0; i < points.size(); ++i) {
            graph[i].resize(points.size());
            for (size_t j = 0; j < points.size(); ++j) {
                if (points[i].distance(points[j]) < d){
                    graph[i][j] = 1;
                }
                else graph[i][j] = 0;
            }
    }
    set<int> used;

    for (size_t i = 0; i < points.size(); ++i) {
        if (used.count(i) != 0) continue;
        Cluster cluster;
        stack<int> stack;
        stack.push(i);
        while (!stack.empty()) {
            int cur = stack.top();
            stack.pop();
            for (size_t j = 0; j < points.size(); ++j) {
                if (graph[cur][j] && used.count(j) == 0) {
                    used.insert(j);
                    cluster.add_point(points[j]);
                    stack.push(j);
                } 
            }
        }
        clusters.push_back(cluster);
    }
    return clusters;
}

vector<Cluster> ClusterFinder::find_clusters_with_k_means_algorithm(vector<Point> points) {
    vector<Cluster> clusters;
    vector<Cluster> res_clusters;
    double mu = 0;
    double min_mu = 1e9;
    double min_k = 0;
    for (size_t k = 1; k < points.size(); ++k) {
        clusters.resize(k);
        mu = k_means(k, clusters, points);
        if (mu < min_mu) {
            min_mu = mu;
            min_k = k;
            for (size_t i = 0; i < clusters.size(); ++i) {
                vector<Point> p;
                for (size_t j = 0; j < clusters[i].points.size(); ++i) {
                    p.push_back(clusters[i].points[j]);
                }
                res_clusters.push_back(Cluster(p));
            }
        }
    }
    return res_clusters;
}

double ClusterFinder::k_means(int k, vector<Cluster>& clusters, vector<Point>& points) {
    vector<int> points_cluster(points.size());
    while (1) {
        for (size_t i = k; i < points.size(); ++i) {
            for (int j = 0; j < k; ++j) {
                return 0;
            }
        }
    }
}

vector<Cluster> ClusterFinder::find_clusters_with_spanning_tree_algorithm(vector<Point> points) {
    vector<Cluster> clusters;/*
    vector<vector<double>> distances_matrix(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        distances_matrix[i].resize(points.size());
        for (size_t j = i; j < points.size(); ++j) {
            distances_matrix[i][j] = -1;
            distances_matrix[j][i] = -1;
        }
    }
    vector<Point> cur_points;
    cur_points.push_back(points[0]);
    set<int> used;
    double min_d;
    int min_i, min_j;

    for (size_t i = 0; i < points.size(); ++i) {
        if (used.count(i) != 0) continue;
        min_d = 1e9;
        for (size_t j = 0; j < cur_points.size(); ++j) {
            for (size_t k = 0; k < points.size(); ++k) {
                if (k != j && k != i && used.count(k) == 0) {
                    if (cur_points[j].distance(points[k]) < min_d) {
                        min_d = cur_points[j].distance(points[k]);
                        min_i = k;
                        min_j = j;
                    }
                }
            }
        }
        cur_points.push_back(points[min_i]);
        used.insert(min_i);
        distances_matrix[min_i][min_j] = min_d;
        distances_matrix[min_j][min_i] = min_d;
        used.insert(i);
    }
    */
    points.size();
    return clusters;
}