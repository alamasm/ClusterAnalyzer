#include "SpanningTreeAlgorithm.h"

void SpanningTreeAlgorithm::find_clusters(vector<Point> points) {
    vector<vector<double>> distances_matrix;
    if (spanning_tree_ == NULL) distances_matrix = spanning_tree(points).first;
    else distances_matrix = spanning_tree_->first;
    for (int i = 0; i < n_clusters; ++i) {
        double max_d = 0;
        int max_j = 0, max_k = 0;
        for (int j = 0; j < distances_matrix.size(); ++j) {
            for (int k = i + 1; k < distances_matrix[j].size(); ++k) {
                if (distances_matrix[j][k] > 0 && distances_matrix[j][k] > max_d) {
                    max_d = distances_matrix[j][k];
                    max_j = j;
                    max_k = k;
                }
            }
        }

        distances_matrix[max_j][max_k] = -1;
        distances_matrix[max_k][max_j] = -1;
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
                if (distances_matrix[cur][j] >= 0 && used.count(j) == 0) {
                    used.insert(j);
                    cluster.add_point(points[j]);
                    stack.push(j);
                } 
            }
        }
        clusters.push_back(cluster);
    }
}

pair<vector<vector<double>>, vector<Point>> SpanningTreeAlgorithm::spanning_tree(vector<Point> points) {
    vector<vector<double>> graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        graph[i].resize(points.size());
        for (size_t j = 0; j < points.size(); ++j) {
            graph[i][j] = -1;
        }
        graph[i][i] = 0;
    }
    vector<int> cur_points;
    cur_points.push_back(0);
    set<int> used;
    used.insert(0);
    double min_d;
    size_t min_i = 0, min_j = 0;

    while (used.size() < points.size()){
        if (used.size()%100 == 0) cout << used.size() << endl;
        min_d = 1e9;
        for (size_t j = 0; j < cur_points.size(); ++j) {
            for (size_t k = 0; k < points.size(); ++k) {
                if (used.count(k) == 0) {
                    if (points[cur_points[j]].distance(points[k]) < min_d) {
                        min_d = points[cur_points[j]].distance(points[k]);
                        min_i = k;
                        min_j = cur_points[j];
                    }
                }
            }
        }
        cur_points.push_back(min_i);
        
        used.insert(min_i);
        graph[min_i][min_j] = points[min_i].distance(points[min_j]);
        graph[min_j][min_i] = points[min_i].distance(points[min_j]);
    }
    return make_pair(graph, points);
}