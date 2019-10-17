#include <vector>
#include <stack>
#include <set>
#ifndef cluster_include
#define cluster_include
#include "Cluster.h"
#endif
#ifndef plane_include
#define plane_include
#include "Plane.h"
#endif

class ClusterFinder
{
    public:
    static const int ALGORITHM_WAVE = 0;
    static const int ALGORITHM_SPANNING_TREE = 1;
    static const int ALGORITHM_K_MEANS = 2;
    ClusterFinder(int algorithm, int d = 0);
    vector<Cluster> find_clusters(Plane plane);
    private:
    int d;
    int algorithm;
    vector<Cluster> find_clusters_with_wave_algorithm(vector<Point> points);
    vector<Cluster> find_clusters_with_k_means_algorithm(vector<Point> points);
    double k_means(int k, vector<Cluster>& clusters, vector<Point>& points);
    vector<Cluster> find_clusters_with_spanning_tree_algorithm(vector<Point> points);
};