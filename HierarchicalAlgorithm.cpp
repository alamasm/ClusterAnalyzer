#include "HierarchicalAlgorithm.h"

void HierarchicalAlgorithm::find_clusters(vector<Point> points) {
    vector<vector<double>> graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        graph[i].resize(points.size());
        for (size_t j = 0; j < points.size(); ++j) {
            graph[i][j] = points[i].distance(points[j]);
        }
    }

    for (int i = 0; i < points.size(); ++i) {
        clusters.push_back(Cluster({points[i]}));
    }

    for (int k_clusters = 0; k_clusters < points.size() - k - 1; k_clusters++) {
        double min_d = 1e9;
        int min_i = 0, min_j = 0;
        for (int i = 0; i < graph.size(); ++i) {
            for (int j = i + 1; j < graph[i].size(); ++j) {
                if (graph[i][j] < min_d) {
                    min_d = graph[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }
        //cout << min_i << " " << min_j << endl;
        union_clusters(min_i, min_j);
        recalc_graph(graph, min_i, min_j);
        clusters.erase(clusters.begin() + min_j);
        //cout <<"sizes: " << clusters.size() << " " << graph.size() << endl;
    }
}

void HierarchicalAlgorithm::union_clusters(int i, int j) {
    //cout << "before: " << clusters[i].points.size() << " " << clusters[j].points.size() << endl;
    for (int k = 0; k < clusters[j].points.size(); ++k) {
        clusters[i].add_point(clusters[j].points[k]);
    }
    //cout << "after: " << clusters[i].points.size() << " " << clusters[j].points.size() << endl;
    //clusters.erase(clusters.begin() + j);
}

void HierarchicalAlgorithm::recalc_graph(vector<vector<double>>& graph, int i, int j) {
    for (int k = 0; k < graph.size(); ++k) {
        if (k == j) continue;
        double d = min(graph[j][k], graph[i][k]);
        graph[i][k] = d;
        graph[k][i] = d;
    }
    graph.erase(graph.begin() + j);
    for (int i = 0; i < graph.size(); ++i) {
        graph[i].erase(graph[i].begin() + j);
    }
}

double HierarchicalAlgorithm::distance_1(vector<Point> a, vector<Point> b) {
    double min = 1e9;
    for (auto p1 : a) {
        for (auto p2 : b) {
            if (p1.distance(p2) < min) {
                min = p1.distance(p2);
            }
        }
    }
    return min;
}