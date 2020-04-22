#include "FORELAlgorithm.h"

void FORELAlgorithm::find_clusters(vector<Point> points) {
    Point point0 = points[0];
    set<int> used;
    bool going = 1;
    vector<int> cur_points;
    while (used.size() < points.size() - 1) {
        while (going) {
            cur_points.clear();
            cur_points.resize(0);
            for (int i = 1; i < points.size(); ++i) {
                if (used.count(i) == 0 && points[i].distance(point0) < R) cur_points.push_back(i);
            }
            Point new_point0 = Point(0, 0);
            for (int i = 0; i < cur_points.size(); ++i) {
                new_point0 += points[cur_points[i]];
            }
            
            new_point0 += point0;
            new_point0 /= (cur_points.size() + 1);
            going = 0;
            if (!(new_point0 == point0)) going = 1;
            point0 = new_point0;
        }

        Cluster cluster;
        for (int i = 0; i < cur_points.size(); ++i) {
            used.insert(cur_points[i]);
            cluster.add_point(points[cur_points[i]]);
        }
        clusters.push_back(cluster);
        going = 1;

        for (int i = 0; i < points.size(); ++i) {
            if (used.count(i) == 0 && !(points[i] == point0)) {
                point0 = points[i];
                break;
            }
        }
    }
}
