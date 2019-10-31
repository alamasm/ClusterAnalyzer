#include "ClusterFinder.h"

#include <iostream>
#include <random>
#define N_ITER 3

ClusterFinder::ClusterFinder(int algorithm, int d) {
    this->algorithm = algorithm;
    this->d = d;
}

vector<Cluster> ClusterFinder::find_clusters(Plane plane) {
    switch (algorithm) {
        case ALGORITHM_WAVE:
            return find_clusters_with_wave_algorithm(plane.get_points());
        case ALGORITHM_SPANNING_TREE:
            return find_clusters_with_spanning_tree_algorithm(plane.get_points());
        case ALGORITHM_K_MEANS:
            return find_clusters_with_k_means_algorithm(plane.get_points());
    }
    return vector<Cluster>();
}

vector<Cluster> ClusterFinder::find_clusters_with_wave_algorithm(vector<Point> points) {
    vector<vector<int>> graph(points.size());
    vector<Cluster> clusters;
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
    return clusters;
}

vector<Cluster> ClusterFinder::find_clusters_with_k_means_algorithm(vector<Point> points) {
    vector<Cluster> clusters;
    vector<Cluster> cur_clusters;
    vector<Cluster> res_clusters;
    double mu = 0;
    double min_mu = 1e9;
    double min_k = 0;
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

        cout << k << " " << min_mu_cur << endl;
        if (min_mu_cur < min_mu) {
            min_mu = min_mu_cur;
            min_k = k;
            copy_clusters(clusters, res_clusters);
        }
    }
    cout << min_k << endl;
    return res_clusters;
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
            cluster_centers[points_cluster[i]] += points[i];
        }
        
        for (long unsigned i = 0; i < k; ++i){
             cluster_centers[i] /= clusters_size[i];
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

vector<Cluster> ClusterFinder::find_clusters_with_spanning_tree_algorithm(vector<Point> points) {
    vector<Cluster> clusters;
    vector<vector<double>> distances_matrix = spanning_tree(points).first;
    
    return clusters;
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
    size_t min_i = 0, min_j = 0, last_i = 0;

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
        last_i = min_i;
    }
    return make_pair(graph, points);
}
