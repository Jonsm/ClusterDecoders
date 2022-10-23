//
//  clustering_helper_naive.cpp
//  ClusterDecoders
//
//  Created by Jon on 1/2/22.
//
#include <iostream>
#include <math.h>
#include "clustering_helper_naive.hpp"

using namespace std;
clustering_helper_naive::clustering_helper_naive(int w, int h) :
w(w),
h(h),
syndromes(w*h),
syndromes_to_index(boost::extents[w][h]),
visited(w*h),
marked(w*h)
{
    fill(syndromes_to_index.data(), syndromes_to_index.data()+syndromes_to_index.num_elements(), -1);
}

//add a new syndrome
void clustering_helper_naive::flip(int x, int y) {
    if (syndromes_to_index[x][y] == -1) {
        syndromes_to_index[x][y] = n_syndromes;
        syndromes[n_syndromes] = pair<int,int>(x,y);
        n_syndromes++;
    } else {
        int ind = syndromes_to_index[x][y];
        syndromes_to_index[x][y] = -1;
        if (ind != n_syndromes - 1) {
            syndromes[ind] = syndromes[n_syndromes - 1];
            syndromes_to_index[syndromes[ind].first][syndromes[ind].second] = ind;
        }
        n_syndromes--;
    }
}

//find clusters (return as a pair of x,y) of separation 2^level
void clustering_helper_naive::get_clusters_helper(vector<vector<pair<int,int>>>& clusters, int level) {
    fill(visited.begin(), visited.begin()+n_syndromes, false);
    int r = pow(2, level);
    int n_visited = 0;
    
    while (n_visited < n_syndromes - n_marked) {
        clusters.push_back(vector<pair<int,int>>());
        vector<pair<int,int>>& cluster = clusters.back();
        
        int first_unvisited = 0;
        while (visited[first_unvisited] == true || marked[first_unvisited] == true) {
            first_unvisited++;
        }
        to_visit.push(first_unvisited);
        
        while(!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();
            
            if (!visited[current] && !marked[current]) {
                for (int i = 0; i < n_syndromes; i++) {
                    if (!visited[i] && !marked[i]) {
                        int dist_x = min(abs(syndromes[current].first - syndromes[i].first), w - abs(syndromes[current].first - syndromes[i].first));
                        int dist_y = min(abs(syndromes[current].second - syndromes[i].second), h - abs(syndromes[current].second - syndromes[i].second));
                        if (max(dist_x, dist_y) <= r) {
                            to_visit.push(i);
                        }
                    }
                }
                
                cluster.push_back(syndromes[current]);
                n_visited++;
                visited[current] = true;
            }
        }
    }
}

//sorting method for pairs which sorts by second element (ie y coord)
bool clustering_helper_naive::sort_second(const pair<int,int> &a, const pair<int,int> &b) {
    return (a.second < b.second);
}

//gets the rect for a cluster
void clustering_helper_naive::get_size(vector<std::pair<int, int>> &cluster, vector<int> &rect) {
    sort(cluster.begin(), cluster.end());
    int max_x_gap = (cluster[0].first - cluster.back().first + w) % w - 1;
    int tl_x = cluster[0].first;
    if (cluster[0].first == cluster.back().first) {
        max_x_gap = w - 1;
    } else {
        for (int i = 1; i < cluster.size(); i++) {
            int x_gap = cluster[i].first - cluster[i-1].first - 1;
            if (x_gap > max_x_gap) {
                max_x_gap = x_gap;
                tl_x = cluster[i].first;
            }
        }
    }
    
    sort(cluster.begin(), cluster.end(), sort_second);
    int max_y_gap = (cluster[0].second - cluster.back().second + h) % h - 1;
    int tl_y = cluster[0].second;
    if (cluster[0].second == cluster.back().second) {
        max_y_gap = h - 1;
    } else {
        for (int i = 1; i < cluster.size(); i++) {
            int y_gap = cluster[i].second - cluster[i-1].second - 1;
            if (y_gap > max_y_gap) {
                max_y_gap = y_gap;
                tl_y = cluster[i].second;
            }
        }
    }
    
    int rect_w = w - max_x_gap;
    int rect_h = h - max_y_gap;
    rect[0] = tl_x;
    rect[1] = tl_y;
    rect[2] = rect_w;
    rect[3] = rect_h;
}

//get all rects
void clustering_helper_naive::get_rects_helper(vector<vector<pair<int,int>>>& clusters, vector<vector<int>>& rects) {
    for (vector<pair<int,int>>& cluster : clusters) {
        rects.push_back(vector<int>(4));
        get_size(cluster, rects.back());
    }
}

//find clusters (return as a pair of x,y) of separation 2^level, and rects bounding each cluster
void clustering_helper_naive::get_clusters(vector<vector<pair<int,int>>>& clusters, vector<vector<int>>& rects, int level) {
    get_clusters_helper(clusters, level);
    get_rects_helper(clusters, rects);
}

int clustering_helper_naive::get_n_syndromes() {
    return n_syndromes;
}

//mark a group of syndromes as belonging to cluster, so we don't need to cluster them again
void clustering_helper_naive::mark_clustered(int x, int y) {
    int ind = syndromes_to_index[x][y];
    marked[ind] = true;
    n_marked++;
}

void clustering_helper_naive::unmark_all() {
    fill(marked.begin(), marked.begin() + n_syndromes, false);
    n_marked = 0;
}

//reset between decoding runs
void clustering_helper_naive::clear() {
    fill(marked.begin(), marked.end(), false);
    n_marked = 0;
    fill(syndromes_to_index.data(), syndromes_to_index.data()+syndromes_to_index.num_elements(), -1);
    n_syndromes = 0;
    fill(visited.begin(), visited.end(), false);
}
