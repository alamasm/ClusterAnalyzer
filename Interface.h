#include <iostream>
#include <string>
#ifndef control_include
#define control_include
#include "Control.h"
#endif

class Interface
{
    public:
    Interface(Control control);
    void start();
    private:
    Control control;
    ifstream input;
    ofstream out;
    ofstream log_out;
    int points_saved = 0;
    int from_file = 0;
    void parse();
    void print_distances(vector<vector<double>> g, ofstream& out);
    void print_points(ofstream &out, vector<Point> points, int group);
    void print_clusters(ofstream& out, ofstream& out_eighen, vector<Cluster> clusters);
    void print_result(ofstream &out);
    void print_spanning_tree(ofstream &out, pair<vector<vector<double>>, vector<Point>> g);
    void log(string s);
    vector<Group> get_groups_from_file(ifstream &in);
};
