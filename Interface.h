#include <iostream>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <string.h>
#include <string>
#include <queue>
#ifndef control_include
#define control_include
#include "Control.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/poll.h>
#include <netinet/in.h>

#endif

class Interface
{
    public:
    Interface(Control control);
    void start();
    void parse(string s);
    void log(string s);
    private:
    Control control;
    ifstream input;
    ofstream out;
    ofstream log_out;
    int points_saved = 0;
    int from_file = 0;
    
    void print_distances(vector<vector<double>> g, ofstream& out);
    void print_points(ofstream &out, vector<Point> points, int group);
    void print_clusters(ofstream& out, ofstream& out_eighen, vector<Cluster> clusters);
    void print_result(ofstream &out);
    void print_spanning_tree(ofstream &out, pair<vector<vector<double>>, vector<Point>> g);
    void print_em_ellipses(ofstream &out_gnu, ofstream &out_em_data, string filename_em_data);
    void print_em_animation(ofstream& out_gnu, ofstream& out_em_animation_data, string filename_em_animation_data);


    

    vector<Group> get_groups_from_file(ifstream &in);
};
