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
        /*if (s == "find_clusters") {
            int algorithm;
            int d;
            *in >> algorithm;
            *in >> d;
            log("finding clusters with " + to_string(algorithm) + " algorithm and d = " + to_string(d));
            control.find_clusters(algorithm, d);
        }*/
        if (s == "find_clusters_wave") {
            double d;
            *in >> d;
            log("finding clusters with wave algorithm and d =" + to_string(d));
            control.find_clusters_wave(d);
        }
        if (s == "find_clusters_k_means") {
            int k;
            *in >> k;
            log("finding clusters with k-means algorithm and k = " + to_string(k));
            control.find_clusters_k_means(k);
        }
        if (s == "find_clusters_spanning_tree") {
            int k;
            *in >> k;
            log("finding clusters with spanning tree algorithm with k = " + to_string(k));
            control.find_clusters_spanning_tree(k);
        }
        if (s == "find_clusters_hierarchical") {
            int k;
            *in >> k;
            log("finding clusters with hierarchical algorithm and k = " + to_string(k));
            control.find_clusters_hierarchical(k);
        }
        if (s == "find_clusters_forel") {
            double R;
            *in >> R;
            log("finding clusters with FOREL algorithm and R = " + to_string(R));
            control.find_clusters_forel(R);
        }
        if (s == "find_clusters_dbscan") {
            int min_pts;
            double eps;
            *in >> min_pts >> eps;
            log("finding clusters with DBSCAN algorithm and min_pts, eps = " + to_string(min_pts) + ", " + to_string(eps));
            control.find_clusters_dbscan(min_pts, eps);
        }
        if (s == "find_clusters_em") {
            int k;
            *in >> k;
            log("finding clusters with EM-Algorithm and k = " + to_string(k));
            control.find_clusters_em(k);
        }
        if (s == "print_clusters") {
            string filename;
            *in >> filename;
            int n;
            *in >> n;
            ofstream out;
            ofstream out_gnu;
            ofstream out_eighen;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            string filename_eighen = filename.substr(0, filename.find_last_of('.')) + "_eighen.txt";
            out_gnu.open(filename_gnu);
            out_gnu << "set palette model RGB defined (0 \"red\",1 \"blue\", 2 \"green\", 3 \"yellow\", 4 \"orange\", 5 \"black\", 6 \"violet\")\n";
            out_gnu << "plot '"+ filename + "' using 1:2:3 notitle with points pt 2 palette, " + "'" + filename_eighen + "' using 1:2:3:4 with vectors filled head lw 3";
            out.open("out/" + filename);
            out_eighen.open("out/" + filename_eighen);
            log("printing clusters to " + filename);
            print_clusters(out, out_eighen, control.get_clusters(n));
            out.close();
            out_gnu.close();
            out_eighen.close();
        }
        if (s == "print_clusters_em") {
            string filename;
            *in >> filename;
            ofstream out;
            ofstream out_gnu;
            ofstream out_animation_gnu;
            ofstream out_eighen;
            ofstream out_em_data;
            ofstream out_em_animation_data;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            string filename_animation_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_animation_gnu.plt";
            string filename_eighen = filename.substr(0, filename.find_last_of('.')) + "_eighen.txt";
            string filename_em_data = filename.substr(0, filename.find_last_of('.')) + "_em_data.txt";
            string filename_em_animation_data = filename.substr(0, filename.find_last_of('.')) + "_animation_data.txt";
            string filename_gif = filename.substr(0, filename.find_last_of('.')) + "_animation.gif";
            out_gnu.open(filename_gnu);
            out_eighen.open("out/" + filename_eighen);
            out_em_data.open("out/" + filename_em_data);
            out_em_animation_data.open("out/" + filename_em_animation_data);
            out_animation_gnu.open(filename_animation_gnu);

            int k = control.get_em_data().first;
            int s = control.get_em_animation_data().second.size();
            out_gnu << "set palette model RGB defined (0 \"red\",1 \"blue\", 2 \"green\", 3 \"yellow\", 4 \"orange\", 5 \"black\", 6 \"violet\")\n";    
            out_gnu << "plot '"+ filename + "' using 1:2:3 notitle with points pt 2 palette, " + "'" + filename_eighen + "' using 1:2:3:4 with vectors filled head lw 3, " + "'" + filename_em_data + "' using 1:2:3:4:5 with ellipses";
            out_animation_gnu << "set term gif animate optimize delay 100 size 600, 600 background \"#ffeedf\" crop" << endl;
            out_animation_gnu << "set output '" << filename_gif << "'" << endl;
            out_animation_gnu << "set size square" << endl;
            out_animation_gnu << "do for [i=0:" << s << "] {" << endl;
            out_animation_gnu << "set palette model RGB defined (0 \"red\",1 \"blue\", 2 \"green\", 3 \"yellow\", 4 \"orange\", 5 \"black\", 6 \"violet\")\n";    
            out_animation_gnu << "plot '"+ filename + "' using 1:2:3 notitle with points pt 2 palette, " + "'" + filename_eighen + "' using 1:2:3:4 with vectors filled head lw 3, " + "'" + filename_em_animation_data + "' every::i*" << k << "::(i+1)*" << k << "-1 using 1:2:3:4:5 with ellipses";
            out_animation_gnu << "}" << endl;
            out.open("out/" + filename);
            log("printing em clusters to " + filename);
            print_clusters(out, out_eighen, control.get_clusters(control.em_index));
            print_em_ellipses(out_gnu, out_em_data, filename_em_data);
            print_em_animation(out_animation_gnu, out_em_animation_data, filename_em_animation_data);
            
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
        if (s == "print_spanning_tree_distances_histagram") {
            string filename;
            *in >> filename;
            ofstream out;
            ofstream out_gnu;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            out_gnu.open(filename_gnu);
            out.open("out/" + filename);
            vector<vector<double>> g = control.get_spanning_tree().first; 
            print_distances(g, out);

            out_gnu << ("binwidth=0.1\n"
                       "bin(x,width)=width*floor(x/width) + width/2.0\n"
                       "set boxwidth binwidth\n"
                       "plot '"+filename+"' using (bin($1,binwidth)):(1.0) smooth freq with boxes");

            log("printing spanning tree histagram to out/" + filename);
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
    for (size_t i = 0; i < g.first.size(); ++i) {
        for (size_t j = i; j < g.first[i].size(); ++j) {
            if (g.first[i][j] >= 0) out << g.second[j].x << " " << g.second[j].y << endl;
        }
        out << endl;
    }
}

void Interface::print_em_ellipses(ofstream& out_gnu, ofstream& out_em_data, string filename_em_data) {
    pair<int, EM_step_data> data = control.get_em_data();
    for (int i = 0; i < data.first; ++i) {
        //out_em_data << data.second.second[i].x << " " << data.second.second[i].y << " " << data.second.first[i].first.x << " " << data.second.first[i].first.y << " "
         //<< data.second.first[i].second << " " << i << endl;
        out_em_data << data.second.centers[i].x << " " << data.second.centers[i].y << " " << data.second.diams[i].x << " " << data.second.diams[i].y << " " << data.second.angles[i] << " " << i << endl;
    }
}

void Interface::print_em_animation(ofstream& out_gnu, ofstream& out_em_animation_data, string filename_em_animation_data) {
    pair<int, vector<EM_step_data>> data = control.get_em_animation_data();
    for (int i = 0; i < data.second.size(); ++i) {
        for (int j = 0; j < data.first; ++j) {
            out_em_animation_data << data.second[i].centers[j].x << " " << data.second[i].centers[j].y << " " << data.second[i].diams[j].x << " " << data.second[i].diams[j].y << " " << data.second[i].angles[j] << " " << j << endl;
        }
    }
}

void Interface::print_distances(vector<vector<double>> g, ofstream& out) {
    for (size_t i = 0; i < g.size(); ++i) {
        for (size_t j = i + 1; j < g.size(); ++j) {
            if (g[i][j] > 0) out << g[i][j] << endl;
        }
    }
}

void Interface::print_points(ofstream& out, vector<Point> points, int group) {
    for (auto p : points) {
        out << p.x << " " << p.y << " " << group << endl;
    }
}

void Interface::print_clusters(ofstream& out, ofstream& out_eighen, vector<Cluster> clusters) {
    for (size_t i = 0; i < clusters.size(); ++i) {
        print_points(out, clusters[i].points, i);
        auto eighen = clusters[i].get_eighen_vectors();
        out_eighen << eighen.first.x << " " << eighen.first.y << " " << eighen.second.first.x << " " << eighen.second.first.y << endl;
        out_eighen << eighen.first.x << " " << eighen.first.y << " " << eighen.second.second.x << " " << eighen.second.second.y << endl;
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
