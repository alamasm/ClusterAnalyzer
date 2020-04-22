#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#endif
class HierarchicalAlgorithm : public Algorithm {
    public:
    int k;
    HierarchicalAlgorithm(int k): k(k) {};
    void find_clusters(vector<Point> points);
    private:
    void union_clusters(int i, int j);
    void recalc_graph(vector<vector<double>>& graph, int i, int j);
    double distance_1(vector<Point> a, vector<Point> b);
};