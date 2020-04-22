#ifndef ClusterFinderDef
#define ClusterFinderDef
#include "ClusterFinder.h"
#endif
using namespace std;
class Algorithm {
    public:
    vector<Cluster> clusters;
    virtual void find_clusters(vector<Point>) = 0;
};