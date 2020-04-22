#ifndef AlgorithmDef
#define AlgorithmDef
#include "Algorithm.h"
#endif
class WaveAlgorithm : public Algorithm {
    public:
    double d;
    WaveAlgorithm(double d) : d(d) {};
    void find_clusters(vector<Point> points);
};