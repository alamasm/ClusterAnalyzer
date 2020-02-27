#include "WaveAlgorithm.h"

void WaveAlgorithm::find_clusters(vector<Point> points) {
    vector<vector<int>> graph(points.size());
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
}
