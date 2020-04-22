#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#include <random>

using namespace std;
#endif
class EMAlgorithm : public Algorithm {
    public:
    vector<vector<vector<long double>>> sigma;
    vector<Point> mu;
    vector<pair<Point, double>> get_eighen();
    int k;
    EMAlgorithm(int k): k(k) {};
    void find_clusters(vector<Point>);
    private:
    bool is_changing(long double llh);
    long double mu_distance(vector<Point> mu);
    long double sigma_distance(vector<vector<vector<long double>>> sigma);
    vector<vector<vector<long double>>> prev_sigma;
    vector<Point> prev_mu;
    long double prev_dist;
    long double prev_llh;
    bool first = 1;
};