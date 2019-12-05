#include "ClusterFinder.h"

#include <iostream>
#include <random>
#include <algorithm>
#define N_ITER 3

ClusterFinder::ClusterFinder(int algorithm, double d) {
    this->algorithm = algorithm;
    this->d = d;
}

void ClusterFinder::find_clusters_with_wave_algorithm(vector<Point> points) {
    vector<vector<int>> graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
            graph[i].resize(points.size());
            for (size_t j = 0; j < points.size(); ++j) {
                if (points[i].distance(points[j]) < d){
                    graph[i][j] = 1;
                }
                else graph[i][j] = 0;
            }
    }
    set<int> used;

    for (size_t i = 0; i < points.size(); ++i) {
        if (used.count(i) != 0) continue;
        Cluster cluster;
        stack<int> stack;
        stack.push(i);
        while (!stack.empty()) {
            int cur = stack.top();
            stack.pop();
            for (size_t j = 0; j < points.size(); ++j) {
                if (graph[cur][j] && used.count(j) == 0) {
                    used.insert(j);
                    cluster.add_point(points[j]);
                    stack.push(j);
                } 
            }
        }
        clusters.push_back(cluster);
    }
}

void ClusterFinder::find_clusters_with_k_means_algorithm(vector<Point> points, int k) {
    vector<Cluster> clusters;
    vector<Cluster> cur_clusters;
    vector<Cluster> res_clusters;
    double mu = 0;
    double min_mu = 1e9;
    double min_k = 0;
    double delta = 0;

    if (k != -1) {
        double min_mu_cur = 1e9;
        for (int i = 0; i < N_ITER; ++i) {
            cur_clusters.clear();
            cur_clusters.resize(k);
            mu = k_means(k, cur_clusters, points);
            if (mu < min_mu_cur) {
                min_mu_cur = mu;
                copy_clusters(cur_clusters, clusters);
            }
        }
        this->clusters = clusters;
        return;
    }

    for (size_t k = 1; k < 7; ++k) {
        clusters.clear();
        clusters.resize(k);

        double min_mu_cur = 1e9;
        for (int i = 0; i < N_ITER; ++i) {
            cur_clusters.clear();
            cur_clusters.resize(k);
            mu = k_means(k, cur_clusters, points);
            if (mu < min_mu_cur) {
                min_mu_cur = mu;
                copy_clusters(cur_clusters, clusters);
            }
        }

        cout << k << " " << min_mu_cur << " " << fabs(min_mu_cur - min_mu) << endl;
        if (fabs(min_mu_cur - min_mu) > delta) {
            min_mu = min_mu_cur;
            min_k = k;
            delta = fabs(min_mu_cur - min_mu);
            copy_clusters(clusters, res_clusters);
        }
    }
    cout << min_k << endl;
    this->clusters = res_clusters;
    //return res_clusters;
}

double ClusterFinder::k_means(long unsigned k, vector<Cluster>& clusters, vector<Point>& points) {
    vector<int> points_cluster(points.size());
    vector<Point> cluster_centers(k);
    vector<int> clusters_size(k);

    double mu = 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    uniform_int_distribution<> dist(0, points.size());
    for (size_t i = 0; i < k; ++i) {
        clusters_size[i] = 1;
        int a = dist(gen);
        cluster_centers[i] = points[a];
    }

    int going = 1;
    while (going) {
        going = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            double min_d = 1e9;
            int min_j = 0;
            for (long unsigned j = 0; j < k; ++j) {
                if (points[i].distance(cluster_centers[j]) < min_d) {
                    min_d = points[i].distance(cluster_centers[j]);
                    min_j = j;
                }
            }
            if (min_j != points_cluster[i]) going = 1;
            points_cluster[i] = min_j;
        }
        

        for (long unsigned i = 0; i < k; ++i) cluster_centers[i] = Point(0, 0);
        for (size_t i = 0; i < points.size(); ++i) {
            clusters_size[points_cluster[i]]++;
            //cout << i << " " << points.size() << endl;
            cluster_centers[points_cluster[i]] += points[i];
        }
        
        for (long unsigned i = 0; i < k; ++i){
             cluster_centers[i] /= (clusters_size[i]);
             clusters_size[i] = 1;
        }
    }
    for (long unsigned i = 0; i < k; ++i) clusters[i] = Cluster();
    for (size_t i = 0; i < points.size(); ++i) clusters[points_cluster[i]].add_point(points[i]);

    for (long unsigned i = 0; i < k; ++i) {
        for (size_t j = 0; j < clusters[i].points.size(); ++j) {
            //mu += clusters[i].points[j].distance(cluster_centers[i]);
            for (size_t k = j + 1; k < clusters[i].points.size(); ++k) {
                mu += clusters[i].points[j].distance(clusters[i].points[k]);
            }
        }
        for (long unsigned j = i + 1; j < k; ++j) {
            mu += cluster_centers[i].distance(cluster_centers[j]);
       }
    }
    /*
    cout << "CENTERS:" << endl;
    for (long unsigned i = 0; i < k; ++i) {
        cout << cluster_centers[i].x << " " << cluster_centers[i].y << endl;
    }*/
    return mu;
}

void ClusterFinder::copy_clusters(vector<Cluster>& from, vector<Cluster>& to) {
    to.clear();
    for (size_t i = 0; i < from.size(); ++i) {
        vector<Point> p;
        for (size_t j = 0; j < from[i].points.size(); ++j) {
            p.push_back(from[i].points[j]);
        }
        to.push_back(Cluster(p));
    }
}

void ClusterFinder::find_clusters_with_spanning_tree_algorithm(vector<Point> points, int n_clusters, pair<vector<vector<double>>, vector<Point>>* spanning_tree_) {
    vector<vector<double>> distances_matrix;
    if (spanning_tree_ == NULL) distances_matrix = spanning_tree(points).first;
    else distances_matrix = spanning_tree_->first;
    for (int i = 0; i < n_clusters; ++i) {
        double max_d = 0;
        int max_j = 0, max_k = 0;
        for (int j = 0; j < distances_matrix.size(); ++j) {
            for (int k = i + 1; k < distances_matrix[j].size(); ++k) {
                if (distances_matrix[j][k] > 0 && distances_matrix[j][k] > max_d) {
                    max_d = distances_matrix[j][k];
                    max_j = j;
                    max_k = k;
                }
            }
        }

        distances_matrix[max_j][max_k] = -1;
        distances_matrix[max_k][max_j] = -1;
    }

    set<int> used;
    for (size_t i = 0; i < points.size(); ++i) {
        if (used.count(i) != 0) continue;
        Cluster cluster;
        stack<int> stack;
        stack.push(i);
        while (!stack.empty()) {
            int cur = stack.top();
            stack.pop();
            for (size_t j = 0; j < points.size(); ++j) {
                if (distances_matrix[cur][j] >= 0 && used.count(j) == 0) {
                    used.insert(j);
                    cluster.add_point(points[j]);
                    stack.push(j);
                } 
            }
        }
        clusters.push_back(cluster);
    }
}

pair<vector<vector<double>>, vector<Point>> ClusterFinder::spanning_tree(vector<Point> points) {
    vector<vector<double>> graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        graph[i].resize(points.size());
        for (size_t j = 0; j < points.size(); ++j) {
            graph[i][j] = -1;
        }
        graph[i][i] = 0;
    }
    vector<int> cur_points;
    cur_points.push_back(0);
    set<int> used;
    used.insert(0);
    double min_d;
    size_t min_i = 0, min_j = 0;

    while (used.size() < points.size()){
        if (used.size()%100 == 0) cout << used.size() << endl;
        min_d = 1e9;
        for (size_t j = 0; j < cur_points.size(); ++j) {
            for (size_t k = 0; k < points.size(); ++k) {
                if (used.count(k) == 0) {
                    if (points[cur_points[j]].distance(points[k]) < min_d) {
                        min_d = points[cur_points[j]].distance(points[k]);
                        min_i = k;
                        min_j = cur_points[j];
                    }
                }
            }
        }
        cur_points.push_back(min_i);
        
        used.insert(min_i);
        graph[min_i][min_j] = points[min_i].distance(points[min_j]);
        graph[min_j][min_i] = points[min_i].distance(points[min_j]);
    }
    return make_pair(graph, points);
}

void ClusterFinder::find_clusters_with_hierarchical_algorithm(vector<Point> points, int k) {
    vector<vector<double>> graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        graph[i].resize(points.size());
        for (size_t j = 0; j < points.size(); ++j) {
            graph[i][j] = points[i].distance(points[j]);
        }
    }

    for (int i = 0; i < points.size(); ++i) {
        clusters.push_back(Cluster({points[i]}));
    }

    for (int k_clusters = 0; k_clusters < points.size() - k - 1; k_clusters++) {
        double min_d = 1e9;
        int min_i = 0, min_j = 0;
        for (int i = 0; i < graph.size(); ++i) {
            for (int j = i + 1; j < graph[i].size(); ++j) {
                if (graph[i][j] < min_d) {
                    min_d = graph[i][j];
                    min_i = i;
                    min_j = j;
                }
            }
        }
        cout << min_i << " " << min_j << endl;
        union_clusters(min_i, min_j);
        recalc_graph(graph, min_i, min_j);
        clusters.erase(clusters.begin() + min_j);
        cout <<"sizes: " << clusters.size() << " " << graph.size() << endl;
    }
}

void ClusterFinder::union_clusters(int i, int j) {
    cout << "before: " << clusters[i].points.size() << " " << clusters[j].points.size() << endl;
    for (int k = 0; k < clusters[j].points.size(); ++k) {
        clusters[i].add_point(clusters[j].points[k]);
    }
    cout << "after: " << clusters[i].points.size() << " " << clusters[j].points.size() << endl;
    //clusters.erase(clusters.begin() + j);
}

void ClusterFinder::recalc_graph(vector<vector<double>>& graph, int i, int j) {
    for (int k = 0; k < graph.size(); ++k) {
        if (k == j) continue;
        double d = min(graph[j][k], graph[i][k]);
        graph[i][k] = d;
        graph[k][i] = d;
    }
    graph.erase(graph.begin() + j);
    for (int i = 0; i < graph.size(); ++i) {
        graph[i].erase(graph[i].begin() + j);
    }
}

double ClusterFinder::distance_1(vector<Point> a, vector<Point> b) {
    double min = 1e9;
    for (auto p1 : a) {
        for (auto p2 : b) {
            if (p1.distance(p2) < min) {
                min = p1.distance(p2);
            }
        }
    }
    return min;
}

void ClusterFinder::find_clusters_with_forel_algorithm(vector<Point> points, double R) {
    Point point0 = points[0];
    set<int> used;
    bool going = 1;
    vector<int> cur_points;
    while (used.size() < points.size() - 1) {
        while (going) {
            cur_points.clear();
            cur_points.resize(0);
            for (int i = 1; i < points.size(); ++i) {
                if (used.count(i) == 0 && points[i].distance(point0) < R) cur_points.push_back(i);
            }
            Point new_point0 = Point(0, 0);
            for (int i = 0; i < cur_points.size(); ++i) {
                new_point0 += points[cur_points[i]];
            }
            
            new_point0 += point0;
            new_point0 /= (cur_points.size() + 1);
            going = 0;
            if (!(new_point0 == point0)) going = 1;
            point0 = new_point0;
        }

        Cluster cluster;
        for (int i = 0; i < cur_points.size(); ++i) {
            used.insert(cur_points[i]);
            cluster.add_point(points[cur_points[i]]);
        }
        clusters.push_back(cluster);
        going = 1;

        for (int i = 0; i < points.size(); ++i) {
            if (used.count(i) == 0 && !(points[i] == point0)) {
                point0 = points[i];
                break;
            }
        }
    }
}

void ClusterFinder::find_clusters_with_dbscan_algorithm(vector<Point> points, int min_pts, double eps) {
    int n_clusters = 0;
    vector<int> label(points.size(), -1);//labels: 0 - noise, c - cluster

    for (size_t i = 0; i < points.size(); ++i) {
        if (label[i] != -1) {
            continue;
        }
        set<int> neighbors = get_neighbors(points, points[i], eps);
        if (neighbors.size() < min_pts) {
            label[i] = 0;
            continue;
        }
        n_clusters++;
        Cluster cluster;    
        label[i] = n_clusters;
        cluster.add_point(points[i]);
        for (size_t j : neighbors) {
            if (j == i) continue;
            
            if (label[j] == 0){
                 label[j] = n_clusters;
                 cluster.add_point(points[j]);
            }
            if (label[j] != -1){continue;}
            label[j] = n_clusters;
            cluster.add_point(points[j]);
            set<int> neighbors2 = get_neighbors(points, points[j], eps);
            if (neighbors2.size() >= min_pts){
                set_union(neighbors.begin(), neighbors.end(), neighbors2.begin(), neighbors2.end(), inserter(neighbors, neighbors.begin()));
            ;}
        }
        clusters.push_back(cluster);
    }
}

set<int> ClusterFinder::get_neighbors(vector<Point> points, Point p, double eps) {
    set<int> neighbors;
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].distance(p) < eps) neighbors.insert(i);
    }
    return neighbors;
}
