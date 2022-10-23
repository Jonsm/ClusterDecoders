//
//  bias_correction_helper.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/16/21.
//

#include "bias_correction_helper.hpp"
#include <iostream>

//DEPRECATED

using namespace std;

typedef boost::multi_array<int, 3>::array_view<3>::type view;
typedef boost::multi_array_types::index_range range;

//w,h = dims, stabilizer (A/B) = which type of stabilizer to correct errors for
bias_correction_helper::bias_correction_helper(int w, int h, char stabilizer) :
w(w),
h(h),
stabilizer(stabilizer),
corrected_syndromes(boost::extents[w][h]),
local_gauss_laws(boost::extents[w][h][w+h])
{
    get_gauss_laws();
}

/*
 * gets a basis of local gauss laws for boxes up to size w,h. Local gauss law for a box is
 * defined as a product of stabilizers such that stabilizers act within the box only as
 * Pauli Z.
 */
void bias_correction_helper::get_gauss_laws() {
    int max_w = h+1;
    vector<int> r1(max_w);
    vector<int> r2(max_w);
    vector<int>* current_row = &r1;
    vector<int>* next_row = &r2;
    
    int x_dir = -1;
    if (stabilizer == 'B') {
        r1[max_w - 1] = 1;
        x_dir = 1;
    } else {
        r1[0] = 1;
    }
    
    for (int y = 0; y < h; y++) {
        next_row->clear();
        for (int x = 0; x < max_w; x++) {
            (*next_row)[x] = (*current_row)[x] ^ (*current_row)[(x+x_dir + max_w) % max_w];
        }
        vector<int>* tmp = current_row;
        current_row = next_row;
        next_row = tmp;
        
        int y_stab = y;
        if (stabilizer == 'B') {
            y_stab = h - y - 1;
        }
        
        for (int offs = -1*h; offs < w; offs++) {
            for (int x = 0; x < max_w; x++) {
                int x_offs = offs + x;
                if (x_offs >= 0 && x_offs < w) {
                    local_gauss_laws[x_offs][y_stab][offs+h] = (*current_row)[x];
                }
            }
        }
    }
}

// check if a cluster with boundary rect is correctable with only pauli Z errors
bool bias_correction_helper::check_cluster(std::vector<std::pair<int, int> > &cluster, std::vector<int> &rect) {
    tl_x = rect[0];
    tl_y = rect[1];
    rect_w = rect[2];
    rect_h = rect[3];
    for (int x = 0; x < rect_w; x++) {
        for (int y = 0; y < rect_h; y++) {
            corrected_syndromes[x][y] = 0;
        }
    }
    
    for (const pair<int,int>& syndrome : cluster) {
        int x_shift = (syndrome.first - tl_x + w) % w;
        int y_shift = (syndrome.second - tl_y + h) % h;
        corrected_syndromes[x_shift][y_shift] = 1;
    }
    
    return gaussian_elimination(cluster);
}

//attempt to correct the cluster with pauli Z operators. If it is not correctable with Z,
//try to correct the most possible syndromes, and then update cluster to contain only the remaining syndromes
void bias_correction_helper::correct_cluster(std::vector<std::pair<int, int> > &cluster, boost::multi_array<int, 2>& correction_Z) {
    update_cluster(cluster);
    apply_correction(correction_Z);
}

//updates vector cluster to contain only the uncorrected syndromes
void bias_correction_helper::update_cluster(vector<pair<int, int> > &cluster) {
    vector<pair<int,int>> new_cluster;
    for (int i = 0; i < cluster.size(); i++) {
        int shift_x = (cluster[i].first - tl_x + w) % w;
        int shift_y = (cluster[i].second - tl_y + h) % h;
        if (corrected_syndromes[shift_x][shift_y] != 1) {
            new_cluster.push_back(cluster[i]);
        }
    }

    cluster = move(new_cluster);
}

//applies a correction of Pauli Z operators. In practice, evolve syndromes using Sierpinski CA rule.
void bias_correction_helper::apply_correction(boost::multi_array<int, 2> &correction_Z) {
    int x_dir = -1;
    int y_dir = -1;
    if (stabilizer == 'B') {
        x_dir = 1; // check
        y_dir = 1;
    }
    
    for (int y = 0; y < rect_h - 1; y++) {
        int shift_y = (y + tl_y + 1) % h;
        int shift_y_local = y;
        if (stabilizer == 'A') {
            shift_y = (tl_y + rect_h - y - 1) % h;
            shift_y_local = rect_h - y - 1;
        }
        
        for (int x = 0; x < rect_w; x++) {
            int shift_x = (x + tl_x) % w;
            if (corrected_syndromes[x][shift_y_local] == 1) {
                correction_Z[shift_x][shift_y] ^= 1;
                corrected_syndromes[x][shift_y_local + y_dir] ^= 1;
                corrected_syndromes[x+x_dir][shift_y_local + y_dir] ^= 1;
            }
        }
    }
}

//helper for determining gauss laws of a box
view bias_correction_helper::get_gauss_law_view(vector<pair<int, int> > &cluster) {
    int gauss_bbox_min_x = 0;
    int gauss_bbox_max_x = rect_w;
    int gauss_bbox_min_y = 0;
    int gauss_bbox_max_y = rect_h;
    int gauss_offs_min = h - rect_h;
    int gauss_offs_max = h + rect_w;
    if (stabilizer == 'B') {
        gauss_bbox_min_x = w-rect_w;
        gauss_bbox_max_x = w;
        gauss_bbox_min_y = h-rect_h;
        gauss_bbox_max_y = h;
        gauss_offs_min = w - rect_w;
        gauss_offs_max = w + rect_h;
    }
    
    return local_gauss_laws[boost::indices[range(gauss_bbox_min_x, gauss_bbox_max_x)][range(gauss_bbox_min_y,gauss_bbox_max_y)][range(gauss_offs_min,gauss_offs_max)]];
}

/*
 * Try to find the largest set of correctable syndromes (greedily, not guaranteed to be largest).
 * Do this by finding largest set that has even support on each local Gauss law. This is done by
 * writing syndrome's support on gauss laws as a vector in Z2 and doing Gaussian elimination to find a
 * sum of syndromes that is 0 mod 2.
 */
bool bias_correction_helper::gaussian_elimination(vector<pair<int,int>>& cluster) {
    view gauss_law_view = get_gauss_law_view(cluster);
    int gauss_law_range = (int)gauss_law_view.shape()[2];

    row_ops.clear();
    current_gl.resize(gauss_law_range);
    pivots.resize(gauss_law_range);
    fill(pivots.begin(), pivots.end(), -1);
    fill(current_gl.begin(), current_gl.end(), 0);
    
    for (const pair<int,int>& syndrome : cluster) {
        int x = (syndrome.first - tl_x + w) % w;
        int y = (syndrome.second - tl_y + h) % h;
        for (int i = 0; i < gauss_law_range; i++) {
            current_gl[i] ^= gauss_law_view[x][y][i];
        }
    }
    
    int original_pivot = 0;
    while (original_pivot != gauss_law_range && current_gl[original_pivot] == 0) {
        original_pivot++;
    }
    if (original_pivot == gauss_law_range) {
        return true;
    }
    for (int i = original_pivot + 1; i < gauss_law_range; i++) {
        if (current_gl[i] == 1) {
            row_ops.push_back(pair<int,int>(original_pivot, i));
        }
    }
    
    for (int i = 0; i < cluster.size(); i++) {
        int x = (cluster[i].first - tl_x + w) % w;
        int y = (cluster[i].second - tl_y + h) % h;
        
        for (int j = 0; j < gauss_law_range; j++) {
            current_gl[j] = gauss_law_view[x][y][j];
        }
        
        for (const pair<int,int>& row_op : row_ops) {
            current_gl[row_op.second] ^= current_gl[row_op.first];
        }
        
        int pivot = 0;
        while (pivot != gauss_law_range && (pivots[pivot] != -1 || pivot == original_pivot || current_gl[pivot] == 0)) {
            pivot++;
        }

        if (pivot == gauss_law_range) {
            if (current_gl[original_pivot] == 1) {
                for (int j = 0; j < gauss_law_range; j++) {
                    if (current_gl[j] == 1 && j != original_pivot) {
                        pair<int,int>& pivot_coord = cluster[pivots[j]];
                        int shift_x = (pivot_coord.first - tl_x + w) % w;
                        int shift_y = (pivot_coord.second - tl_y + h) % h;
                        corrected_syndromes[shift_x][shift_y] = 0;
                    }
                }
                
                int shift_x_final = (cluster[i].first - tl_x + w) % w;
                int shift_y_final = (cluster[i].second - tl_y + h) % h;
                corrected_syndromes[shift_x_final][shift_y_final] = 0;
                
                return false;
            }
        } else {
            pivots[pivot] = i;
            for (int j = 0; j < gauss_law_range; j++) {
                if (current_gl[j] == 1 && j != pivot) {
                    row_ops.push_back(pair<int,int>(pivot, j));
                }
            }
        }
    }
    
    return false;
}
