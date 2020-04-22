#include "EMAlgorithm.h"
#include <tgmath.h>
#include <iostream>
#define STOP_MIN 5

long double det_2x2(vector<vector<long double>> m) {
    return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

long double get_power(Point p, Point mu, vector<vector<long double>> sigma) {
    //return 0.5 * 1/det_2x2(sigma) * ((mu.y - p.y) * (sigma[0][0] * (p.y - mu.y) + (sigma[0][1] + sigma[1][0]) * (p.x - mu.x)) - sigma[1][1] * (p.x - mu.x) * (p.x - mu.x));
    return -0.5 * ((p.x - mu.x) * (p.x - mu.x) / sigma[0][0] + (p.y - mu.y) * (p.y - mu.y) / sigma[1][1]);
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
    /*if (first == 1) {
        prev_mu = mu;
        prev_sigma = sigma;
        first = 0;
        prev_dist = 0;
        return 1;
    }
    double dist = sigma_distance(sigma) + mu_distance(mu);
    cout << "DIST: " << dist << endl;
    if (fabs(dist - prev_dist) <= STOP_MIN) return 0;
    prev_mu = mu;
    prev_sigma = sigma;
    prev_dist = dist;
    return 1;
    */
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
        double lambda1 = 0.5 * (sigma[i][0][0] + sigma[i][1][1] + sqrt((sigma[i][0][0] + sigma[i][1][1]) * (sigma[i][0][0] + sigma[i][1][1]) - 4 * det_2x2(sigma[i])));
        double lambda2 = 0.5 * (sigma[i][0][0] + sigma[i][1][1] - sqrt((sigma[i][0][0] + sigma[i][1][1]) * (sigma[i][0][0] + sigma[i][1][1]) - 4 * det_2x2(sigma[i])));
        double x1 = -lambda1 / sigma[i][1][0];   
        double x2 = -lambda1 / sigma[i][1][0];  
        double cos = 0;

        if (x1 * x1 + 1 >= x2 * x2 + 1) {
            cos = x1 / sqrt(x1 * x1 + 1);
        } else {
            cos = x2 / sqrt(x2 * x2 + 1);
        }
        cout << "X1 X1 " << x1 << " " << x2 << endl;
        cout << "COS  " << cos << endl;
        
        res.push_back(make_pair(Point(2 * max(lambda1, lambda2), 2 * min(lambda1, lambda2)), acos(cos) * 180));
    }
    return res;
}

void EMAlgorithm::find_clusters(vector<Point> points) {
    //init
    vector<long double> w(k);
    for (int i  = 0; i < k; ++i) {
        w[i] = 1.0 / k;
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
            mu[i].x += 1/points.size() * points[l].x;
            mu[i].y += 1/points.size() * points[l].y;
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
    while (is_changing(llh)){
        //E-step
        cout << "E" << endl;
        cout << "BEFORE:" << endl;
        for (int j = 0; j < k; ++j) {
            cout << "j = " << j << endl;
            print_sigma(sigma[j]);
        }
        llh = 0;
        vector<long double> sump(points.size());
        vector<vector<long double>> P_1(points.size());
        //vector<vector<long double>> P_2(k);
        //for (int j = 0; j < k; ++j) P_2[j].resize(points.size());
        vector<Point> new_mu(k);
        vector<long double> new_w(k);
        vector<vector<vector<long double>>> new_sigma(k);
        for (int i = 0; i < k; ++i) {
            new_mu[i].x = 0;
            new_mu[i].y = 0;
            new_w[i] = 0;
        }
        for (int i = 0; i < k; ++i) {
            new_sigma[i].resize(2);
            for (int j = 0; j < 2; ++j) {
                new_sigma[i][j].resize(2);
                for (int l = 0; l < 2; ++l) {
                    new_sigma[i][j][l] = 0;
                }
                sigma[i][j][j] = 1;
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
                new_w[j] += P_1[i][j] / sump[i];
            }
            llh += log(sump[i]);
            //cout << "sump " << sump[i] << endl;
        }
        //M-step
        cout << "M" << endl;
        for (int j = 0; j < k; ++j) {
            mu[j].x = new_mu[j].x / new_w[j];
            mu[j].y = new_mu[j].y / new_w[j];
            //cout << "new_w " << new_w[j] << endl;
            //cout << "mu " << mu[j].x << " " << mu[j].y << endl;
            for (int i = 0; i < points.size(); ++i) {
                new_sigma[j][0][0] += P_1[i][j] / sump[i] * (points[i].x - mu[j].x) * (points[i].x - mu[j].x);
                //new_sigma[j][0][1] += P_1[i][j] / sump[i] * (points[i].x - mu[j].x) * (points[i].y - mu[j].y);
                //new_sigma[j][1][0] += P_1[i][j] / sump[i] * (points[i].x - mu[j].x) * (points[i].y - mu[j].y);
                new_sigma[j][0][1] = 0;
                new_sigma[j][1][0] = 0;
                new_sigma[j][1][1] += P_1[i][j] / sump[i] * (points[i].y - mu[j].y) * (points[i].y - mu[j].y);
            }
           
            sigma[j][0][0] = new_sigma[j][0][0] / points.size();
            sigma[j][0][1] = new_sigma[j][0][1] / points.size();
            sigma[j][1][0] = new_sigma[j][1][0] / points.size();
            sigma[j][1][1] = new_sigma[j][1][1] / points.size();
            cout << "mew sigma j: " << j << endl;
            print_sigma(sigma[j]);
            w[j] = new_w[j] / points.size();
            cout << "new mu j: " << j << " " << mu[j].x << " " << mu[j].y << endl;
        }
        cout << "AFTER:" << endl;
        for (int j = 0; j < k; ++j) {
            cout << "j = " << j << endl;
            print_sigma(sigma[j]);
        }
    }
    Cluster cluster;
    for (int i = 0; i < points.size(); ++i) {
        cluster.add_point(points[i]);
    }
    clusters.push_back(cluster);
}

