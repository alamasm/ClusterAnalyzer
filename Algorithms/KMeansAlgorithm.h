#include "Algorithm.h"
class KMeansAlgorithm : public Algorithm {
    public:
    int k;
    KMeansAlgorithm(int k): k(k) {};
    void find_clusters(vector<Point> points);
    private:
    double k_means(long unsigned k, vector<Cluster>& clusters, vector<Point>& points);
    void copy_clusters(vector<Cluster>& from, vector<Cluster>& to);
};