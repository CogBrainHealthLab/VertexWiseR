#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector TFCE_cpp_impl(NumericVector data, int tail, IntegerMatrix edgelist) {
    int N = data.size();
    int E = edgelist.nrow();
    
    NumericVector tfce_sum(N, 0.0);
    if (N == 0 || E == 0) return tfce_sum;
    
    double max_score = 0.0;
    std::vector<int> signs;
    
    if (tail == 2) {
        signs = {-1, 1};
        for (int i = 0; i < N; ++i) {
            if (!NumericVector::is_na(data[i])) {
                double abs_val = std::abs(data[i]);
                if (abs_val > max_score) max_score = abs_val;
            }
        }
    } else if (tail == 1) {
        signs = {1};
        for (int i = 0; i < N; ++i) {
            if (!NumericVector::is_na(data[i])) {
                if (data[i] > max_score) max_score = data[i];
            }
        }
    } else if (tail == -1) {
        signs = {-1};
        for (int i = 0; i < N; ++i) {
            if (!NumericVector::is_na(data[i])) {
                double neg_val = -data[i];
                if (neg_val > max_score) max_score = neg_val;
            }
        }
    } else {
        Rcpp::stop("Invalid tail parameter. Must be -1, 1, or 2.");
    }
    
    if (max_score <= 0.0) return tfce_sum;
    
    double step = max_score / 100.0;
    std::vector<double> score_threshs(100);
    for (int i = 0; i < 100; ++i) {
        score_threshs[i] = step * (i + 1);
    }
    
    std::vector<std::pair<int, int>> active_edges;
    active_edges.reserve(E);
    for (int i = 0; i < E; ++i) {
        int u = edgelist(i, 0) - 1;
        int v = edgelist(i, 1) - 1;
        if (u >= 0 && u < N && v >= 0 && v < N) {
            active_edges.push_back({u, v});
        }
    }
    
    std::vector<int> parent(N);
    std::vector<int> sz(N);
    std::vector<int> gen(N, 0);
    std::vector<int> visited(N, 0);
    
    std::vector<int> root_to_cid(N, 0);
    std::vector<int> root_to_cid_gen(N, 0);
    std::vector<int> cluster_sizes(N); 
    
    int current_gen = 0;
    
    auto get_parent = [&](int u) {
        int root = u;
        while (parent[root] != root) root = parent[root];
        int curr = u;
        while (curr != root) {
            int next = parent[curr];
            parent[curr] = root;
            curr = next;
        }
        return root;
    };
    
    for (int sign : signs) {
        std::vector<double> val(N);
        for (int i = 0; i < N; ++i) {
            val[i] = NumericVector::is_na(data[i]) ? -1e300 : data[i] * sign;
        }
        
        std::vector<std::pair<int, int>> curr_active_edges = active_edges;
        
        for (int t = 0; t < 100; ++t) {
            double thresh = score_threshs[t];
            double thresh_sq = thresh * thresh;
            
            std::vector<std::pair<int, int>> next_active_edges;
            next_active_edges.reserve(curr_active_edges.size());
            for (auto& e : curr_active_edges) {
                if (val[e.first] >= thresh && val[e.second] >= thresh) {
                    next_active_edges.push_back(e);
                }
            }
            curr_active_edges = std::move(next_active_edges);
            
            if (curr_active_edges.empty()) continue; 
            
            current_gen++;
            int num_clusters = 0;
            
            for (auto& e : curr_active_edges) {
                int u = e.first;
                if (gen[u] != current_gen) { parent[u] = u; sz[u] = 1; gen[u] = current_gen; }
                int v = e.second;
                if (gen[v] != current_gen) { parent[v] = v; sz[v] = 1; gen[v] = current_gen; }
                
                int root_u = get_parent(u);
                int root_v = get_parent(v);
                if (root_u != root_v) {
                    if (sz[root_u] < sz[root_v]) std::swap(root_u, root_v);
                    parent[root_v] = root_u;
                    sz[root_u] += sz[root_v];
                }
            }
            
            for (auto& e : curr_active_edges) {
                int r = get_parent(e.first);
                if (root_to_cid_gen[r] != current_gen) {
                    root_to_cid_gen[r] = current_gen;
                    if (sz[r] > 1) {
                        num_clusters++;
                        root_to_cid[r] = num_clusters;
                        cluster_sizes[num_clusters - 1] = sz[r];
                    } else {
                        root_to_cid[r] = 0;
                    }
                }
            }
            
            if (num_clusters == 0) continue;
            
            for (auto& e : curr_active_edges) {
                int u = e.first;
                if (visited[u] != current_gen) {
                    visited[u] = current_gen;
                    int r = get_parent(u);
                    if (root_to_cid[r] > 0) {
                        tfce_sum[u] += sign * cluster_sizes[root_to_cid[r] - 1] * thresh_sq;
                    }
                }
                
                int v = e.second;
                if (visited[v] != current_gen) {
                    visited[v] = current_gen;
                    int r = get_parent(v);
                    if (root_to_cid[r] > 0) {
                        tfce_sum[v] += sign * cluster_sizes[root_to_cid[r] - 1] * thresh_sq;
                    }
                }
            }
        }
    }
    
    return tfce_sum;
}