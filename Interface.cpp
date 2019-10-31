#include <iostream>
#include "Interface.h"

using namespace std;

Interface::Interface(Control control) : out("out/points.txt"), log_out("out/log.txt"){
    this->control = control;
}

void Interface::start() {
    cout << "Read commands from file 1 from konsole 0:" << endl;
    cin >> from_file;
    if (from_file) {
        string filename;
        cout << "Enter filename:" << endl;
        cin >> filename;
        filename = "in/" + filename;
        input.open(filename);
    }
    if (from_file) {
        cout << "Reading commands from file..." << endl;
    } else {
        cout << "Enter commands here:" << endl;
    }
    parse();
}

void Interface::parse() {
    istream* in;
    string s;
    if (from_file) {
        in = &input;
    } else {
        in = &cin;
    }
    while (1) {
        *in >> s;
        if (s.find_first_of("//") == 0) continue;
        if (s == "end") {
            print_result(out);
            if (!points_saved) {
                ofstream points_save;
                points_save.open("save/points_save.txt");
                print_result(points_save);
                log("saving points to save/points_save.txt");
            }
            log("ending...");
            log_out.close();
            out.close();
            break;
        }
        if (s == "create_group") {
            int n;
            double x0, y0, x_disp, y_disp;
            *in >> n >> x0 >> y0 >> x_disp >> y_disp;
            log("creating group with n = " + to_string(n) + " x0 = " + to_string(x0) + " y0 = " + to_string(y0) + " x_disp = "
             + to_string(x_disp) + " y_disp = " + to_string(y_disp));
            control.create_group(n, x0, y0, x_disp, y_disp);
        }
        if (s == "rotate_center") {
            int n;
            double alpha;
            *in >> n >> alpha;
            log("rotating group " + to_string(n) + " relatively to center on " + to_string(alpha) + " radians");
            control.rotate_group_relatively_to_center(n, alpha);
        }
        if (s == "rotate_origin") {
            int n;
            double alpha;
            *in >> n >> alpha;
            log("rotating group " + to_string(n) + " relatively to origin on " + to_string(alpha) + " radians");
            control.rotate_group_relatively_to_origin(n, alpha);
        }
        if (s == "find_clusters") {
            int algorithm;
            int d;
            *in >> algorithm;
            *in >> d;
            log("finding clusters with " + to_string(algorithm) + " algorithm and d = " + to_string(d));
            control.find_clusters(algorithm, d);
        }
        if (s == "print_clusters") {
            string filename;
            *in >> filename;
            int n;
            *in >> n;
            ofstream out;
            ofstream out_gnu;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            out_gnu.open(filename_gnu);
            out_gnu << "set palette model RGB defined (0 \"red\",1 \"blue\", 2 \"green\", 3 \"yellow\", 4 \"orange\", 5 \"black\", 6 \"violet\")\n";
            out_gnu << "plot '"+ filename + "' using 1:2:3 notitle with points pt 2 palette";
            out.open("out/" + filename);
            log("printing clusters to " + filename);
            print_clusters(out, control.get_clusters(n));
            out.close();
            out_gnu.close();
        }
        if (s == "print_spanning_tree") {
            string filename;
            *in >> filename;
            ofstream out;
            ofstream out_gnu;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            out_gnu.open(filename_gnu);
            pair<vector<vector<double>>, vector<Point>> res = control.get_spanning_tree();
            out_gnu << "plot for [j=0:" + to_string(res.second.size()) + "] '" + filename + "' index j u 1:2 with lines";
            out.open("out/" + filename);
            log("printing spanning tree to " + filename);
            print_spanning_tree(out, res);
            out.close();
            out_gnu.close();
        }
        if (s == "read_points") {
            string filename;
            *in >> filename;
            ifstream in_points;
            if (filename == "-1") {
                in_points.open("save/points_save.txt");
                filename = "points_save.txt";
            } else {
                in_points.open("save/" + filename);
            }
            log("reading points from save/" + filename);
            control.set_groups(get_groups_from_file(in_points));
            in_points.close();
        }
        if (s == "save_points") {
            string filename;
            *in >> filename;
            log("saving points to save/" + filename);
            ofstream points_save;
            points_save.open("save/" + filename);
            print_result(points_save);
            points_save.close();
            points_saved = 1;
        }
    }
}

void Interface::print_spanning_tree(ofstream &out, pair<vector<vector<double>>, vector<Point>> g) {
    for (int i = 0; i < g.first.size(); ++i) {
        for (int j = i; j < g.first[i].size(); ++j) {
            if (g.first[i][j] >= 0) out << g.second[j].x << " " << g.second[j].y << endl;
        }
        out << endl;
    }
}

void Interface::print_points(ofstream& out, vector<Point> points, int group) {
    for (auto p : points) {
        out << p.x << " " << p.y << " " << group << endl;
    }
}

void Interface::print_clusters(ofstream& out, vector<Cluster> clusters) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        print_points(out, clusters[i].points, i);
    }
}

void Interface::print_result(ofstream& out) {
    vector<Group> groups = control.get_groups();

    for (size_t i = 0; i < groups.size(); ++i) {
        print_points(out, groups[i].points, i);
    }
}

void Interface::log(string s) {
    log_out << s << endl;
}

vector<Group> Interface::get_groups_from_file(ifstream &in) {
    int group = 0;
    int cur_group = 0;
    double cur_x = 0, cur_y = 0;
    vector<Group> groups;
    vector<Point> points;
    while (in >> cur_x >> cur_y >> cur_group) {
        if (cur_group != group) {
            groups.push_back(Group(points));
            group = cur_group;
            points.resize(0);
        }
        points.push_back(Point(cur_x, cur_y));
    }
    groups.push_back(Group(points));
    return groups;
}