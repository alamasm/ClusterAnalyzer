#include "Algorithm.h"
class DBScanAlgorithm : public Algorithm {
    public:
    int min_pts;
    double eps;
    DBScanAlgorithm(int min_pts, double eps): min_pts(min_pts), eps(eps){};
    void find_clusters(vector<Point>);
    private:
    set<int> get_neighbors(vector<Point> points, Point p, double eps);
};