#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#endif
class FORELAlgorithm : public Algorithm {
    public:
    double R;
    FORELAlgorithm(double R): R(R) {};
    void find_clusters(vector<Point> points);
};