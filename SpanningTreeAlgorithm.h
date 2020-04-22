#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#endif
class SpanningTreeAlgorithm : public Algorithm {
    public:
    int n_clusters;
    pair<vector<vector<double>>, vector<Point>>* spanning_tree_;
    SpanningTreeAlgorithm(int n_clusters, pair<vector<vector<double>>, vector<Point>>* spanning_tree_) : n_clusters(n_clusters), spanning_tree_(spanning_tree_) {};
    void find_clusters(vector<Point> points);
    pair<vector<vector<double>>, vector<Point>> spanning_tree(vector<Point> points);
};