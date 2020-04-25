#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#include <random>

using namespace std;
#endif
struct EM_step_data {
    vector<Point> centers;
    vector<double> angles;
    vector<Point> diams;
};
class EMAlgorithm : public Algorithm {
    public:
    vector<EM_step_data> animation;
    EM_step_data get_res_data();
    int k;
    EMAlgorithm(int k): k(k) {};
    void find_clusters(vector<Point>);
    private:
    vector<pair<Point, double>> get_eighen();
    vector<vector<vector<long double>>> sigma;
    vector<Point> mu;
    bool is_changing(long double llh);
    long double mu_distance(vector<Point> mu);
    long double sigma_distance(vector<vector<vector<long double>>> sigma);
    vector<vector<vector<long double>>> prev_sigma;
    vector<Point> prev_mu;
    long double prev_dist;
    long double prev_llh;
    bool first = 1;
};