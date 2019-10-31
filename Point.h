class Point
{
    public:
    double x, y;
    Point operator+= (Point b);
    Point operator/= (double x);
    Point(double x, double y);
    //Point(const Point& point);
    Point();
    double distance(Point p);
};
