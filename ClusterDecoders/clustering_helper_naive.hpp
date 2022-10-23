//
//  clustering_helper_naive.hpp
//  ClusterDecoders
//
//  Created by Jon on 1/2/22.
//

/*
 * Helper class for clustering syndromes. Finds clusters of syndromes where all syndromes have at most
 * 2^level distance to the closest other syndrome in the cluster. Works by simply computing pairwise
 * distances for all pairs. This has poor scaling, but in practice is faster than the method proposed
 * by Bravyi and Haah for small temperatures, since the number of syndromes will always be small.
 */

#ifndef clustering_helper_naive_hpp
#define clustering_helper_naive_hpp

#include <stdio.h>
#include <vector>
#include <queue>
#include "boost/multi_array.hpp"

class clustering_helper_naive {
public:
    clustering_helper_naive(int w, int h);
    void flip(int x, int y);
    void get_clusters(std::vector<std::vector<std::pair<int,int>>>& clusters, std::vector<std::vector<int>>& rects, int level);
    int get_n_syndromes();
    void mark_clustered(int x, int y);
    void unmark_all();
    void clear();
    
private:
    const int w;
    const int h;
    std::vector<std::pair<int,int>> syndromes;
    boost::multi_array<int,2> syndromes_to_index;
    int n_syndromes = 0;
    std::queue<int> to_visit;
    std::vector<bool> visited;
    int n_marked = 0;
    std::vector<bool> marked;
    
    void get_clusters_helper(std::vector<std::vector<std::pair<int,int>>>& clusters, int level);
    void get_rects_helper(std::vector<std::vector<std::pair<int,int>>>& clusters, std::vector<std::vector<int>>& rects);
    void get_size(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    
    static bool sort_second(const std::pair<int,int> &a, const std::pair<int,int> &b);
};

#endif /* clustering_helper_naive_hpp */
