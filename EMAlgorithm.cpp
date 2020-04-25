#include "EMAlgorithm.h"
#include <tgmath.h>
#include <iostream>
#define STOP_MIN 0.0001

long double det_2x2(vector<vector<long double>> m) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

long double get_power(Point p, Point mu, vector<vector<long double>> sigma) {
    return (-0.5 / det_2x2(sigma)) * (sigma[0][0] * (p.y - mu.y) * (p.y - mu.y) - (sigma[0][1] + sigma[1][0]) * (p.x - mu.x) * (p.y - mu.y) + sigma[1][1] * (p.x - mu.x) * (p.x - mu.x));
    //return -0.5 * ((p.x - mu.x) * (p.x - mu.x) / sigma[0][0] + (p.y - mu.y) * (p.y - mu.y) / sigma[1][1]);
}

EM_step_data EMAlgorithm::get_res_data() {
    EM_step_data res;
    vector<pair<Point, double>> eighen = get_eighen();
    res.centers = mu;
    for (int i = 0; i < k; ++i) {
        res.diams.push_back(eighen[i].first);
        res.angles.push_back(eighen[i].second);
    }
    return res;
}

long double EMAlgorithm::mu_distance(vector<Point> mu) {
    double res = 0;
    for (int i = 0; i < k; ++i) {
        res += (prev_mu[i].x - mu[i].x) * (prev_mu[i].x - mu[i].x);
        res += (prev_mu[i].y - mu[i].y) * (prev_mu[i].y - mu[i].y);
    }
    //cout << mu[0].x << endl;
    return sqrt(res);
}

void print_sigma(vector<vector<long double>> sigma) {
    cout << sigma[0][0] << " " << sigma[0][1] << endl;
    cout << sigma[1][0] << " " << sigma[1][1] << endl;
}

long double EMAlgorithm::sigma_distance(vector<vector<vector<long double>>> sigma) {
    long double res = 0;
    for (int i = 0; i < k; ++i) {
        res += (prev_sigma[i][0][0] - sigma[i][0][0]) * (prev_sigma[i][0][0] - sigma[i][0][0]);
        res += (prev_sigma[i][0][1] - sigma[i][0][1]) * (prev_sigma[i][0][1] - sigma[i][0][1]);
        res += (prev_sigma[i][1][0] - sigma[i][1][0]) * (prev_sigma[i][1][0] - sigma[i][1][0]);
        res += (prev_sigma[i][1][1] - sigma[i][1][1]) * (prev_sigma[i][1][1] - sigma[i][1][1]);
        res /= (sigma[i][0][0] * sigma[i][0][0] + sigma[i][1][0] * sigma[i][1][0] + sigma[i][0][1] * sigma[i][0][1] + sigma[i][1][1] * sigma[i][1][1]);
        //res /= det_2x2(prev_sigma[i]);
    }
    //cout << sigma[0][0][0] << endl;
    return res;
}

bool EMAlgorithm::is_changing(long double llh) {
   cout << "DIST: " << fabs(llh - prev_llh) << endl;
   if (first == 1) {
       prev_llh = llh;
       first = 0;
       return 1;
   }
   if (fabs(llh - prev_llh) <= STOP_MIN) return 0;
   prev_llh = llh;
   return 1;
}

vector<pair<Point, double>> EMAlgorithm::get_eighen() {
    vector<pair<Point ,double>> res;
    for (int i = 0; i < k; ++i){ 
        long double r = sigma[i][0][1] / sqrt(sigma[i][0][0] * sigma[i][1][1]);
        long double phi;
        cout << "R: " << r << endl;
        if (fabs(sigma[i][0][0] - sigma[i][1][1]) > 0.00001) {
            phi = 0.5 * atan(2 * r * sigma[i][0][0] * sigma[i][1][1] / (sigma[i][0][0] * sigma[i][0][0] - sigma[i][1][1] * sigma[i][1][1]));
        } else phi = 0;
        cout << phi << endl;

        //res.push_back(make_pair(Point(2 * sqrt(max(sigma[i][0][0], sigma[i][1][1])), 2 * sqrt(min(sigma[i][0][0], sigma[i][1][1]))),180 - phi / M_PI * 180));
        res.push_back(make_pair(Point(4 * sqrt(sigma[i][0][0]), 4 * sqrt(sigma[i][1][1])) , 180 - phi / M_PI * 180));
    }
    return res;
}

void EMAlgorithm::find_clusters(vector<Point> points) {
    //init
    vector<long double> w(k);
    vector<long double> N(k);
    for (int i  = 0; i < k; ++i) {
        w[i] = 1.0 / k;
        N[i] = points.size() / k;
    }
    sigma = vector<vector<vector<long double>>>(k);
    //vector<vector<vector<double>>> sigma(k);
    for (int i = 0; i < k; ++i) {
        sigma[i].resize(2);
        for (int j = 0; j < 2; ++j) {
            sigma[i][j].resize(2);
            for (int l = 0; l < 2; ++l) {
                sigma[i][j][l] = 0;
            }
            sigma[i][j][j] = 1;
        }
    } 
    mu = vector<Point>(k); 
    for (int i = 0; i < k; ++i) {
        mu[i].x = 0;
        mu[i].y = 0;
        
        for (int l = 0; l < points.size(); ++l) {
            mu[i].x += (1.0 / points.size()) * points[l].x;
            mu[i].y += (1.0 / points.size()) * points[l].y;
        }
    }
    double max_d_x = 0, max_d_y = 0;
    for (int i = 0; i < points.size(); ++i) {
        for (int j = i + 1; j < points.size(); ++j) {
            if (points[i].x - points[j].x > max_d_x) max_d_x = points[i].x - points[j].x;
            if (points[i].y - points[j].y > max_d_y) max_d_y = points[i].y - points[j].y;
        }
    }
    for (int i = 0; i < k; ++i) {
        std::random_device rd;
        std::mt19937 gen(rd());;
        uniform_real_distribution<double> unif_x(-max_d_x / 2, max_d_x / 2);
        uniform_real_distribution<double> unif_y(-max_d_y / 2, max_d_y / 2);
        mu[i].x += unif_x(gen);
        mu[i].y += unif_y(gen);
        cout << "MU " << mu[i].x << " " << mu[i].y << endl;
    }
    long double llh = 0;
    vector<long double> sump(points.size());
        vector<vector<long double>> P_1(points.size());
        //vector<vector<long double>> P_2(k);
        //for (int j = 0; j < k; ++j) P_2[j].resize(points.size());
        vector<Point> new_mu(k);
        vector<long double> new_w(k);
        vector<long double> new_N(k);
        vector<vector<vector<long double>>> new_sigma(k);
    while (is_changing(llh)){
        //E-step
        animation.push_back(get_res_data());
        cout << "E" << endl;
        cout << "BEFORE:" << endl;
        for (int j = 0; j < k; ++j) {
            cout << "j = " << j << endl;
            print_sigma(sigma[j]);
        }
        llh = 0;
        
        for (int i = 0; i < k; ++i) {
            new_mu[i].x = 0;
            new_mu[i].y = 0;
            new_N[i] = 0;
            new_w[i] = 0;
        }
        for (int i = 0; i < k; ++i) {
            new_sigma[i].resize(2);
            for (int j = 0; j < 2; ++j) {
                new_sigma[i][j].resize(2);
                for (int l = 0; l < 2; ++l) {
                    new_sigma[i][j][l] = 0;
                }
                //new_sigma[i][j][j] = 1;
            }
            
        } 

        for (int i = 0; i < points.size(); ++i) {
            sump[i] = 0;
            P_1[i].resize(k);
            for (int j = 0; j < k; ++j) {
                P_1[i][j] = (w[j] / (2 * M_PI * sqrt(det_2x2(sigma[j])))) * expl(get_power(points[i], mu[j], sigma[j]));
                //cout << "power " << get_power(points[i], mu[j], sigma[j]) << endl;
                //cout << "P_1 " << P_1[i][j] << endl;
                //cout << "w " << w[j] << endl;
                //cout << "den " << sqrt(2 * M_PI * det_2x2(sigma[j])) << endl;
                sump[i] += P_1[i][j];
            }
            for (int j = 0; j < k; ++j) {
                new_mu[j].x += (P_1[i][j] / sump[i]) * points[i].x;
                new_mu[j].y += (P_1[i][j] / sump[i]) * points[i].y;
                //new_w[j] += (P_1[i][j] / sump[i]);
                new_N[j] += (P_1[i][j] / sump[i]);
            }
            llh += log(sump[i]);
            //cout << "sump " << sump[i] << endl;
        }
        //M-step
        cout << "M" << endl;
        for (int j = 0; j < k; ++j) {
            mu[j].x = new_mu[j].x / (new_N[j]);
            mu[j].y = new_mu[j].y / (new_N[j]);
            //mu[j].x = new_mu[j].x / 
            //cout << "new_w " << new_w[j] << endl;
            //cout << "mu " << mu[j].x << " " << mu[j].y << endl;
            for (int i = 0; i < points.size(); ++i) {
                new_sigma[j][0][0] += (1 / new_N[j]) * (P_1[i][j] / sump[i]) * (points[i].x - mu[j].x) * (points[i].x - mu[j].x);
                new_sigma[j][0][1] += (1 / new_N[j]) * (P_1[i][j] / sump[i]) * (points[i].x - mu[j].x) * (points[i].y - mu[j].y);
                new_sigma[j][1][0] += (1 / new_N[j]) * (P_1[i][j] / sump[i]) * (points[i].x - mu[j].x) * (points[i].y - mu[j].y);
                new_sigma[j][1][1] += (1 / new_N[j]) * (P_1[i][j] / sump[i]) * (points[i].y - mu[j].y) * (points[i].y - mu[j].y);
            }
            sigma[j][0][0] = new_sigma[j][0][0];
            sigma[j][0][1] = new_sigma[j][0][1];
            sigma[j][1][0] = new_sigma[j][1][0];
            sigma[j][1][1] = new_sigma[j][1][1];
            cout << "mew sigma j: " << j << endl;
            print_sigma(sigma[j]);
            w[j] = new_N[j] / points.size();
            N[j] = new_N[j];
            cout << "new mu j: " << j << " " << mu[j].x << " " << mu[j].y << endl;
        }
        cout << "AFTER:" << endl;
        for (int j = 0; j < k; ++j) {
            cout << "j = " << j << endl;
            print_sigma(sigma[j]);
        }
    }
    animation.push_back(get_res_data());
    for (int i = 0; i < k; ++i) {
        clusters.push_back(Cluster());
    }
    for (int i = 0; i < points.size(); ++i) {
        long double max_p = 0;
        int max_p_j = 0;
        for (int j = 0; j < k; ++j) {
            if (P_1[i][j] > max_p) {
                max_p = P_1[i][j];
                max_p_j = j;
            }
        }
        clusters[max_p_j].add_point(points[i]);
    }
    /*
    for (int i = 0; i < points.size(); ++i) {
        cluster.add_point(points[i]);
    }
    clusters.push_back(cluster);
    */
}