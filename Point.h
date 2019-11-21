class Point
{
    public:
    const double eps = 0.0001;
    double x, y;
    Point operator+= (Point b);
    Point operator/= (double x);
    bool operator == (Point b);
    Point operator= (Point b);
    Point(double x, double y);
    //Point(const Point& point);
    Point();
    double distance(Point p);
};
