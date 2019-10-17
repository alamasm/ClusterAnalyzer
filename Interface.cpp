#include <iostream>
#include "Interface.h"

Interface::Interface(Control control) : out("out/points.txt"){
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
        if (s == "end") {
            control.print_result(out);
            break;
        }
        if (s == "create_group") {
            int n;
            double x0, y0, x_disp, y_disp;
            *in >> n >> x0 >> y0 >> x_disp >> y_disp;
            control.create_group(n, x0, y0, x_disp, y_disp);
        }
        if (s == "rotate_center") {
            int n;
            double alpha;
            *in >> n >> alpha;
            control.rotate_group_relatively_to_center(n, alpha);
        }
        if (s == "rotate_origin") {
            int n;
            double alpha;
            *in >> n >> alpha;
            control.rotate_group_relatively_to_origin(n, alpha);
        }
        if (s == "find_clusters") {
            int algorithm;
            int d;
            *in >> algorithm;
            *in >> d;
            control.find_clusters(algorithm, d);
        }
        if (s == "print_clusters") {
            string filename;
            *in >> filename;
            filename = "out/" + filename;
            ofstream out;
            ofstream out_gnu;
            string filename_gnu = "out/" + filename.substr(0, filename.find_last_of('.')) + "_gnu.plt";
            out_gnu.open(filename_gnu);
            out_gnu << "set palette model RGB defined (0 \"red\",1 \"blue\", 2 \"green\", 3 \"yellow\", 4 \"orange\", 5 \"black\", 6 \"violet\")\n";
            out_gnu << "plot '"+ filename + "' using 1:2:3 notitle with points pt 2 palette";
            out.open(filename);
            control.print_clusters(out);
            out.close();
            out_gnu.close();
        }
    }
}