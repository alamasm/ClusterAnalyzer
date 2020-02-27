#include "Algorithm.h"
class FORELAlgorithm : public Algorithm {
    public:
    double R;
    FORELAlgorithm(double R): R(R) {};
    void find_clusters(vector<Point> points);
};