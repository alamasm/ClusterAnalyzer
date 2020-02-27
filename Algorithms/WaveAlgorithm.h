#include "Algorithm.h"
class WaveAlgorithm : public Algorithm {
    public:
    double d;
    WaveAlgorithm(double d) : d(d) {};
    void find_clusters(vector<Point> points);
};