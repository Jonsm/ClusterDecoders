//
//  clustering_helper.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

//helper class to find clusters in an array of stabilizer syndromes. Works as described in
//sec VI of https://arxiv.org/pdf/1112.3252.pdf

#ifndef clustering_helper_hpp
#define clustering_helper_hpp

#include <vector>
#include <queue>
#include "boost/multi_array.hpp"

struct extremes {
    int bottom;
    int top;
    int left;
    int right;
};

enum dirs {
    N = 0, NW, W, SW, S, SE, E, NE
};

class clustering_helper {
public:
    clustering_helper(int w, int h);
    
    //find clusters. Outputs them to clusters, rects, where each vector in clusters is a list of
    //coordinates of syndromes in the cluster, and each vector in rects is the bounding box of the
    //cluster. Each vector corresponds to a rect as follows: {top left x, top left y, width, height}.
    void get_clusters(boost::multi_array<int, 2>& syndromes, std::vector<std::vector<std::pair<int,int>>>& clusters, std::vector<std::vector<int>>& rects, int level);
private:
    int w;
    int h;
    int w_coarse;
    int h_coarse;
    int r;
    boost::multi_array<int,2> already_clustered;
    std::queue<std::pair<int,int>> to_visit;
    boost::multi_array<extremes, 2> box_extremes;
    boost::multi_array<int,2> ll_bounds;
    boost::multi_array<int,2> lr_bounds;
    boost::multi_array<int,3> delta;
    std::vector<std::pair<int,int>> dir_pairs;
    
    void get_delta(boost::multi_array<int, 2>& syndromes);
    void get_delta_vh(int x, int y);
    void get_delta_diagonal(boost::multi_array<int, 2>& syndromes, int x, int y);
    void preprocess_box_delta(boost::multi_array<int, 2>& syndromes, int idx, int idy);
    void coord_to_rect(int x, int y, int& top, int& left, int& rect_w, int& rect_h);
    void visit(boost::multi_array<int, 2>& syndromes, std::pair<int,int> coord, std::vector<std::pair<int,int>>& current_cluster);
    void visit_r1(boost::multi_array<int, 2>& syndromes, std::pair<int,int> coord, std::vector<std::pair<int,int>>& current_cluster);
    void visit_r_larger(boost::multi_array<int, 2>& syndromes, std::pair<int,int> coord, std::vector<std::pair<int,int>>& current_cluster);
    void get_size(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
};

#endif /* clustering_helper_hpp */
