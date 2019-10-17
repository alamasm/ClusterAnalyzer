#include <vector>
#ifndef group_include
#define group_include
#include "Group.h"
#endif
class Plane
{
    public:
    vector<Group> groups;
    void add_group(int n, double x0, double y0, double x_dispersion, double y_dispersion);
    int rotate_group_relatively_to_center(size_t group_n, double alpha);
    int rotate_group_relatively_to_origin(size_t group_n, double alpha);
    vector<Point> get_points();
};