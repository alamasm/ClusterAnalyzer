#include "KMeansAlgorithm.h"
#define N_ITER 10
void KMeansAlgorithm::find_clusters(vector<Point> points) {
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

        //cout << k << " " << min_mu_cur << " " << fabs(min_mu_cur - min_mu) << endl;
        if (fabs(min_mu_cur - min_mu) > delta) {
            min_mu = min_mu_cur;
            min_k = k;
            delta = fabs(min_mu_cur - min_mu);
            copy_clusters(clusters, res_clusters);
        }
    }
    //cout << min_k << endl;
    this->clusters = res_clusters;
    //return res_clusters;
}

double KMeansAlgorithm::k_means(long unsigned k, vector<Cluster>& clusters, vector<Point>& points) {
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

void KMeansAlgorithm::copy_clusters(vector<Cluster>& from, vector<Cluster>& to) {
    to.clear();
    for (size_t i = 0; i < from.size(); ++i) {
        vector<Point> p;
        for (size_t j = 0; j < from[i].points.size(); ++j) {
            p.push_back(from[i].points[j]);
        }
        to.push_back(Cluster(p));
    }
}