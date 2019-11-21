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
    static const int ALGORITHM_HIERARCHICAL = 3;

    ClusterFinder(int algorithm, double d = 0);
    pair<vector<vector<double>>, vector<Point>> spanning_tree(vector<Point> points);
    vector<Cluster> clusters;
    //vector<Cluster> find_clusters(Plane plane);
    
    void find_clusters_with_wave_algorithm(vector<Point> points);
    void find_clusters_with_k_means_algorithm(vector<Point> points, int k);
    void find_clusters_with_spanning_tree_algorithm(vector<Point> points, int n_clusters, pair<vector<vector<double>>, vector<Point>> *spanning_tree = NULL);
    void find_clusters_with_hierarchical_algorithm(vector<Point> points, int k);
    void find_clusters_with_forel_algorithm(vector<Point> points, double R);

    static double distance_1(vector<Point> a, vector<Point> b);

    private:
    double d;
    int algorithm;
    
    double k_means(long unsigned k, vector<Cluster>& clusters, vector<Point>& points);
    
    void copy_clusters(vector<Cluster>& from, vector<Cluster>& to);
    void recalc_graph(vector<vector<double>>& graph, int i, int j);
    void union_clusters(int i, int j);
};
