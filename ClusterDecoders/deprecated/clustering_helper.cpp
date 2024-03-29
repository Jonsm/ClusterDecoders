//
//  clustering_helper.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

//DEPRECATED. Works but is not as fast in practive as the naive O(errors^2) method.

#include <iostream>
#include <math.h>
#include "boost/multi_array.hpp"
#include "clustering_helper.hpp"

using namespace std;

typedef boost::multi_array<int, 2>::array_view<2>::type view;
typedef boost::multi_array_types::index_range range;

clustering_helper::clustering_helper(int w, int h) :
w(w),
h(h),
already_clustered(boost::extents[w][h]),
partial_boxes(boost::extents[w][h]),
empty(boost::extents[w][h]),
box_extremes(boost::extents[w][h]),
ll_bounds(boost::extents[w][h]),
lr_bounds(boost::extents[w][h]),
delta(boost::extents[w][h][8])
{
    dir_pairs = vector<pair<int,int>>{pair<int,int>(0,-1),pair<int,int>(-1,-1),pair<int,int>(-1,0),pair<int,int>(-1,1),pair<int,int>(0,1),pair<int,int>(1,1),pair<int,int>(1,0),pair<int,int>(1,-1)};
}

//gets delta(i,j) for regions i,j as described in https://arxiv.org/pdf/1112.3252.pdf. I.e.
//whether they are part of the same cluster.
void clustering_helper::get_delta(boost::multi_array<int, 2> &syndromes) {
    fill(lr_bounds.data(), lr_bounds.data()+lr_bounds.num_elements(), 0);
    fill(ll_bounds.data(), ll_bounds.data()+ll_bounds.num_elements(), 0);
    fill(delta.data(), delta.data()+delta.num_elements(), 0);
    
    for (int x = 0; x < w_coarse; x++) {
        for (int y = 0; y < h_coarse; y++) {
            preprocess_box_delta(syndromes, x, y);
        }
    }
    
    for (int x = 0; x < w_coarse; x++) {
        for (int y = 0; y < h_coarse; y++) {
            get_delta_vh(x,y);
            get_delta_diagonal(syndromes,x,y);
        }
    }
}

void clustering_helper::get_delta_vh_edge(int x, int y) {
    int top_c, left_c, w_c, h_c;
    coord_to_rect(x, y, top_c, left_c, w_c, h_c);
    
    if (partial_boxes[x][y] == E || partial_boxes[x][y] == SE) {
        int x_left = x-1;
        int x_right = (x + 1) % w_coarse;
        int top_r, left_r, w_r, h_r, top_l, left_l, w_l, h_l;
        coord_to_rect(x_left, y, top_l, left_l, w_l, h_l);
        coord_to_rect(x_right, y, top_r, left_r, w_r, h_r);
        
        int x_diff = w_l - box_extremes[x_left][y].right + box_extremes[x_right][y].left + w_c;
        if (x_diff <= r) {
            delta[x][y][W] = 1;
            delta[x][y][E] = 1;
            delta[x_left][y][E] = 1;
            delta[x_right][y][W] = 1;
        }
    }
    
    if (partial_boxes[x][y] == S || partial_boxes[x][y] == SE) {
        int y_top = y - 1;
        int y_bot = (y + 1) % h_coarse;
        int top_t, left_t, w_t, h_t, top_b, left_b, w_b, h_b;
        coord_to_rect(x, y_top, top_t, left_t, w_t, h_t);
        coord_to_rect(x, y_bot, top_b, left_b, w_b, h_b);
        
        int y_diff = h_t - box_extremes[x][y_top].bottom + box_extremes[x][y_bot].top + h_c;
        if (y_diff <= r) {
            delta[x][y][N]=1;
            delta[x][y][S]=1;
            delta[x][y_top][S]=1;
            delta[x][y_bot][N]=1;
        }
    }
}

//checks if regions with same x or y are in the same cluster
void clustering_helper::get_delta_vh(int x, int y) {
    if (empty[x][y]) {
        if (partial_boxes[x][y] != 0) {
            get_delta_vh_edge(x, y);
        }
        return;
    }
    int x_right = (x + 1) % w_coarse;
    int y_bottom = (y + 1) % h_coarse;
    
    int top_c, left_c, w_c, h_c;
    coord_to_rect(x, y, top_c, left_c, w_c, h_c);
    
    int x_diff = w_c - box_extremes[x][y].right + box_extremes[x_right][y].left;
    int y_diff = h_c - box_extremes[x][y].bottom + box_extremes[x][y_bottom].top;
    
    if (!empty[x_right][y] && x_diff <= r) {
        delta[x][y][E] = 1;
        delta[x_right][y][W] = 1;
    }
    if (!empty[x][y_bottom] && y_diff <= r) {
        delta[x][y][S] = 1;
        delta[x][y_bottom][N] = 1;
    }
}

bool clustering_helper::get_delta_diagonal_tr_helper(boost::multi_array<int, 2> &syndromes, int x_bl, int y_bl, int x_tr, int y_tr) {
    if (empty[x_bl][y_bl] || empty[x_tr][y_tr]) {
        return false;
    }
    
    int top_tr,left_tr,w_tr,h_tr,top_bl,left_bl,w_bl,h_bl,ofs_x,ofs_y;
    coord_to_rect(x_tr, y_tr, top_tr, left_tr, w_tr, h_tr);
    coord_to_rect(x_bl, y_bl, top_bl, left_bl, w_bl, h_bl);
    view syndromes_bl = syndromes[boost::indices[range(left_bl,left_bl+w_bl)][range(top_bl,top_bl+h_bl)]];
    view bounds_tr = ll_bounds[boost::indices[range(left_tr,left_tr+w_tr)][range(top_tr,top_tr+h_tr)]];
    
    int x_bl_right = left_bl + w_bl;
    ofs_x = (left_tr - x_bl_right + w) % w;
    int y_tr_bot = top_tr + h_tr;
    ofs_y = (top_bl - y_tr_bot + w) % w;
    
    for (int x = 0; x < w_bl; x++) {
        for (int y = 0; y < h_bl; y++) {
            if (syndromes_bl[x][y] == 1) {
                int tr_x_furthest = max(0, min(r-(w_bl-x)-ofs_x,w_tr-1));
                int tr_y_furthest = min(h_tr-1, max(0, h_tr-(r-y)+ofs_y));
                if (bounds_tr[tr_x_furthest][tr_y_furthest] == 1) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

bool clustering_helper::get_delta_diagonal_tl_helper(boost::multi_array<int, 2> &syndromes, int x_br, int y_br, int x_tl, int y_tl) {
    if (empty[x_br][y_br] || empty[x_tl][y_tl]) {
        return false;
    }
    
    int top_tl,left_tl,w_tl,h_tl,top_br,left_br,w_br,h_br,ofs_x,ofs_y;
    coord_to_rect(x_tl, y_tl, top_tl, left_tl, w_tl, h_tl);
    coord_to_rect(x_br, y_br, top_br, left_br, w_br, h_br);
    view syndromes_br = syndromes[boost::indices[range(left_br,left_br+w_br)][range(top_br,top_br+h_br)]];
    view bounds_tl = lr_bounds[boost::indices[range(left_tl,left_tl+w_tl)][range(top_tl,top_tl+h_tl)]];
    
    int x_tl_right = left_tl+w_tl;
    ofs_x = (left_br - x_tl_right + w) % w;
    int y_tl_bot = top_tl+h_tl;
    ofs_y = (top_br - y_tl_bot + h) % h;
    
    for (int x = 0; x < w_br; x++) {
        for (int y = 0; y < h_br; y++) {
            if (syndromes_br[x][y] == 1) {
                int tl_x_furthest = max(0, w_tl-(r-x)-ofs_x);
                int tl_y_furthest = min(h_tl-1, max(0, h_tl-(r-y)+ofs_y));
                if (bounds_tl[tl_x_furthest][tl_y_furthest] == 1) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

void clustering_helper::get_delta_diagonal_edge(boost::multi_array<int, 2> &syndromes, int x, int y) {
    int x_tr, y_tr, x_tl, y_tl, x_br, y_br, x_bl, y_bl;
    if (partial_boxes[x][y] == S || partial_boxes[x][y] == SE) {
        x_tl = x;
        x_br = (x + 1) % w_coarse;
        y_tl = (y - 1 + h_coarse) % h_coarse;
        y_br = (y + 1) % h_coarse;
        if (get_delta_diagonal_tl_helper(syndromes, x_br, y_br, x_tl, y_tl)) {
            delta[x][y][N] = 1;
            delta[x][y][SE] = 1;
            delta[x_tl][y_tl][S] = 1;
            delta[x_br][y_br][NW] = 1;
        }
        
        x_tr = x;
        x_bl = (x - 1 + w_coarse) % w_coarse;
        y_tr = (y - 1 + h_coarse) % h_coarse;
        y_bl = (y + 1) % h_coarse;
        if (get_delta_diagonal_tl_helper(syndromes, x_bl, y_bl, x_tr, y_tr)) {
            delta[x][y][N] = 1;
            delta[x][y][SW] = 1;
            delta[x_tr][y_tr][S] = 1;
            delta[x_bl][y_bl][NE] = 1;
        }
    }
    
    if (partial_boxes[x][y] == E || partial_boxes[x][y] == SE) {
        x_tl = (x - 1 + w_coarse) % w_coarse;
        x_br = (x + 1) % w_coarse;
        y_tl = (y - 1 + h_coarse) % h_coarse;
        y_br = y;
        if (get_delta_diagonal_tl_helper(syndromes, x_br, y_br, x_tl, y_tl)) {
            delta[x][y][NW] = 1;
            delta[x][y][E] = 1;
            delta[x_tl][y_tl][SE] = 1;
            delta[x_br][y_br][W] = 1;
        }
        
        x_tr = (x + 1) % w_coarse;
        x_bl = (x - 1 + w_coarse) % w_coarse;
        y_tr = (y - 1 + h_coarse) % h_coarse;
        y_bl = y;
        if (get_delta_diagonal_tl_helper(syndromes, x_bl, y_bl, x_tr, y_tr)) {
            delta[x][y][NE] = 1;
            delta[x][y][W] = 1;
            delta[x_tr][y_tr][SW] = 1;
            delta[x_bl][y_bl][E] = 1;
        }
    }
    
    if (partial_boxes[x][y] == SE) {
        x_tl = (x - 1 + w_coarse) % w_coarse;
        x_br = (x + 1) % w_coarse;
        y_tl = (y - 1 + h_coarse) % h_coarse;
        y_br = (y + 1) % h_coarse;
        if (get_delta_diagonal_tl_helper(syndromes, x_br, y_br, x_tl, y_tl)) {
            delta[x][y][NW] = 1;
            delta[x][y][SE] = 1;
            delta[x_tl][y_tl][SE] = 1;
            delta[x_br][y_br][NW] = 1;
        }
        
        x_tr = (x + 1) % w_coarse;
        x_bl = (x - 1 + w_coarse) % w_coarse;
        y_tr = (y - 1 + h_coarse) % h_coarse;
        y_bl = (y + 1) % h_coarse;
        if (get_delta_diagonal_tl_helper(syndromes, x_bl, y_bl, x_tr, y_tr)) {
            delta[x][y][NE] = 1;
            delta[x][y][SW] = 1;
            delta[x_tr][y_tr][SW] = 1;
            delta[x_bl][y_bl][NE] = 1;
        }
    }
}

//checks if regions touching diagonally are in the same cluster
void clustering_helper::get_delta_diagonal(boost::multi_array<int, 2> &syndromes, int idx, int idy) {
    if (empty[idx][idy] && partial_boxes[idx][idy] != 0) {
        get_delta_diagonal_edge(syndromes, idx, idy);
    }
    
    int x_tr = (idx + 1) % w_coarse;
    int x_tl = (idx - 1 + w_coarse) % w_coarse;
    int y_top = (idy - 1 + h_coarse) % h_coarse;
    
    int top_c, left_c, w_c, h_c, top_tr, left_tr, w_tr, h_tr, top_tl, left_tl, w_tl, h_tl;
    coord_to_rect(idx, idy, top_c, left_c, w_c, h_c);
    coord_to_rect(x_tr, y_top, top_tr, left_tr, w_tr, h_tr);
    coord_to_rect(x_tl, y_top, top_tl, left_tl, w_tl, h_tl);
    view syndromes_c = syndromes[boost::indices[range(left_c,left_c+w_c)][range(top_c,top_c+h_c)]];
    view bounds_tr = ll_bounds[boost::indices[range(left_tr,left_tr+w_tr)][range(top_tr,top_tr+h_tr)]];
    view bounds_tl = lr_bounds[boost::indices[range(left_tl,left_tl+w_tl)][range(top_tl,top_tl+h_tl)]];
    
    bool empty_tl = empty[x_tl][y_top];
    bool empty_tr = empty[x_tr][y_top];
    if (empty_tl && empty_tr) {
        return;
    }

    for (int x = 0; x < w_c; x++) {
        for (int y = 0; y < h_c; y++) {
            if (syndromes_c[x][y] == 1) { 
                if (!empty_tl) {
                    int tl_x_furthest = max(0, w_tl-(r-x));
                    int tl_y_furthest = max(0, h_tl-(r-y));
                    if (bounds_tl[tl_x_furthest][tl_y_furthest] == 1) {
                        empty_tl = true;
                        delta[idx][idy][NW] = 1;
                        delta[x_tl][y_top][SE] = 1;
                    }
                }
                if (!empty_tr) {
                    int tr_x_furthest = min(r-(w_c-x),w_tr-1);
                    int tr_y_furthest = max(0, h_tr-(r-y));
                    if (bounds_tr[tr_x_furthest][tr_y_furthest] == 1) {
                        empty_tr = true;
                        delta[idx][idy][NE] = 1;
                        delta[x_tr][y_top][SW] = 1;
                    }
                }
            }
        }
    }
}

//helper for getting delta, checks the region (box), finds location of syndromes in the region and
//neighboring regions
void clustering_helper::preprocess_box_delta(boost::multi_array<int, 2> &syndromes, int idx, int idy) {
    int top, left, box_w, box_h;
    coord_to_rect(idx, idy, top, left, box_w, box_h);
//    cout << "BOX:" << top << "," << left << "," << box_w << "," << box_h << endl;
    view box_view = syndromes[boost::indices[range(left,left+box_w)][range(top,top+box_h)]];

    extremes& extremes = box_extremes[idx][idy];
    extremes.top = box_h - 1;
    extremes.bottom = 0;
    extremes.left = box_w - 1;
    extremes.right = 0;
    bool is_empty = true;
    
    for (int x = 0; x < box_w; x++) {
        for (int y = 0; y < box_h; y++) {
            if (box_view[x][y] == 1) {
                is_empty = false;
                if (x > extremes.right) {
                    extremes.right = x;
                }
                if (x < extremes.left) {
                    extremes.left = x;
                }
                if (y > extremes.bottom) {
                    extremes.bottom = y;
                }
                if (y < extremes.top) {
                    extremes.top = y;
                }
            }
        }
    }
    
    int ll_top = 0;
    int x = 0;
    while (x < box_w && ll_top < box_h) {
        for (int y = box_h - 1; y >= ll_top; y--) {
            if (box_view[x][y] == 1) {
                for (int x2 = x; x2 < box_w; x2++) {
                    for (int y2 = y; y2 >= ll_top; y2--) {
                        ll_bounds[left+x2][top+y2] = 1;
                    }
                }
                ll_top = y + 1;
                    
                break;
            }
        }
        
        x++;
    }
    
    int lr_top = 0;
    x = box_w-1;
    while (x >= 0 && lr_top < box_h) {
        for (int y = box_h-1; y >= lr_top; y--) {
            if (box_view[x][y] == 1) {
                for (int x2 = x; x2 >= 0; x2--) {
                    for (int y2 = y; y2 >= lr_top; y2--) {
                        lr_bounds[left+x2][top+y2] = 1;
                    }
                }
                break;
            }
        }
        
        x--;
    }
    
    if (is_empty) {
        empty[idx][idy] = 1;
    }
}

void clustering_helper::preprocess_partial() {
    for (int x = 0; x < w_coarse; x++) {
        for (int y = 0; y < h_coarse; y++) {
            if (x == w_coarse - 1 && y == h_coarse - 1) {
                partial_boxes[x][y] = SE;
            } else if (x == w_coarse - 1) {
                partial_boxes[x][y] = E;
            } else if (y == h_coarse - 1) {
                partial_boxes[x][y] = S;
            } else {
                partial_boxes[x][y] = 0;
            }
        }
    }
}

/*
finds clusters as in https://arxiv.org/pdf/1112.3252.pdf, i.e. maximal regions of syndromes
with distance no more than 2^level from some other syndrome.
 ** syndromes = input syndromes
 ** clusters = vector to output clusters, where each vector of pairs is a cluster
 ** rects = vector to output rects, where each rect is a bounding box of cluster at same position
    in clusters. rect is a vector of {top left, top right, width, height}
 ** level = max distance of elements from all other elements in cluster
 */
void clustering_helper::get_clusters(boost::multi_array<int, 2> &syndromes, vector<vector<pair<int, int>>> &clusters, vector<vector<int>> &rects, int level) {
    
    r = pow(2, level);
    w_coarse = ceil((float)w / r);
    h_coarse = ceil((float)h / r);
    
    for (int x = 0; x < w_coarse; x++) {
        for (int y = 0; y < h_coarse; y++) {
            already_clustered[x][y] = 0;
            empty[x][y] = 0;
        }
    }
    
    if (r != 1) {
        preprocess_partial();
        get_delta(syndromes);
    }
    
    clusters.push_back(vector<pair<int,int>>());
    for (int x = 0; x < w_coarse; x++) {
        for (int y = 0; y < h_coarse; y++) {
            if (!already_clustered[x][y]) {
                to_visit.push(pair<int,int>(x,y));
                
                while (!to_visit.empty()) {
                    visit(syndromes, to_visit.front(), clusters.back());
                    to_visit.pop();
                }
                
                if (!clusters.back().empty()) {
                    rects.push_back(vector<int>(4));
                    get_size(clusters.back(), rects.back());
                    
                    clusters.push_back(vector<pair<int,int>>());
                }
            }
        }
    }
    clusters.pop_back();
}

//visit a region in BFS
void clustering_helper::visit(boost::multi_array<int, 2> &syndromes, pair<int,int> coord, vector<pair<int,int>>& current_cluster) {
    if (already_clustered[coord.first][coord.second]) {
        return;
    }
    already_clustered[coord.first][coord.second] = 1;
    
    if (r == 1) {
        visit_r1(syndromes, coord, current_cluster);
    } else {
        visit_r_larger(syndromes, coord, current_cluster);
    }
}

//special case of visit for r = 1 (reduces to just checking adjacency)
void clustering_helper::visit_r1(boost::multi_array<int, 2> &syndromes, pair<int, int> coord, vector<pair<int, int> > &current_cluster) {
    if (syndromes[coord.first][coord.second] == 0) {
        return;
    }
    
    current_cluster.push_back(coord);
    
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            int shift_x = (x + coord.first + w) % w;
            int shift_y = (y + coord.second + h) % h;
            
            if ((x != 0 || y != 0) && syndromes[shift_x][shift_y] == 1) {
                to_visit.push(pair<int,int>(shift_x, shift_y));
            }
        }
    }
}

//visit regions larger than 1 site
void clustering_helper::visit_r_larger(boost::multi_array<int, 2> &syndromes, std::pair<int, int> coord, std::vector<std::pair<int, int> > &current_cluster) {
    int top, left, box_w, box_h;
    coord_to_rect(coord.first, coord.second, top, left, box_w, box_h);
    
    if (!empty[coord.first][coord.second]) {
        for (int x = left; x < left + box_w; x++) {
            for (int y = top; y < top + box_h; y++) {
                if (syndromes[x][y] == 1) {
                    current_cluster.push_back(pair<int,int>(x,y));
                }
            }
        }
    }
    
    for (int i = 0; i < 8; i++) {
        if (delta[coord.first][coord.second][i] == 1) {
            int shift_x = (coord.first + dir_pairs[i].first + w_coarse) % w_coarse;
            int shift_y = (coord.second + dir_pairs[i].second + h_coarse) % h_coarse;
            to_visit.push(pair<int,int>(shift_x, shift_y));
        }
    }
}

//gets the rect of region given by coord x,y
void clustering_helper::coord_to_rect(int x, int y, int &top, int &left, int &rect_w, int &rect_h) {
    top = y * r;
    left = x * r;
    int bottom_edge = min(h, top + r);
    int right_edge = min(w, left + r);
    rect_h = bottom_edge - top;
    rect_w = right_edge - left;
}

//sorting method for pairs which sorts by second element (ie y coord)
bool sort_second(const pair<int,int> &a, const pair<int,int> &b) {
    return (a.second < b.second);
}

//gets the rect for a cluster
void clustering_helper::get_size(vector<std::pair<int, int>> &cluster, vector<int> &rect) {
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
