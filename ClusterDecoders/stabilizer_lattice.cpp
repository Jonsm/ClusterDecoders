//
//  StabilizerLattice.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

#include <iostream>
#include "math.h"
#include "boost/multi_array.hpp"
#include "stabilizer_lattice.hpp"

using namespace std;

/*
 **w,h = dims
 **bias = array of {pI,px,py,pz} (currently unused)
 **stabilizer = which type of stabilizer errors to correct, A or B
 **give_up = whether the RG decoder should report failure if clusters are larger than max(w,h)/2
 **reduce_weight = whether decoder should reduce weight of correction by multiplying with stabilizers
   (good for debugging)
 */
stabilizer_lattice::stabilizer_lattice(int w, int h, vector<float> bias, char stabilizer, bool give_up, bool reduce_weight) :
w(w),
h(h),
stabilizer(stabilizer),
errors_Z(boost::extents[w][h]),
errors_not_Z(boost::extents[w][h]),
syndromes(boost::extents[w][h]),
correction_Z(boost::extents[w][h]),
correction_not_Z(boost::extents[w][h]),
clustering_helper(w,h),
reduce_weight(reduce_weight),
give_up(give_up),
bias_correction(w, h, stabilizer)
{
    this->bias = bias;
    
    vector<pair<int,int>> v1 = {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    vector<pair<int,int>> v2 = {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(1,0)};
    if (stabilizer == 'A') {
        offsets_Z = v1;
        offsets_not_Z = v2;
    } else {
        offsets_Z = v2;
        offsets_not_Z = v1;
    }
}

//flips a qubit at site x,y. z = whether Pauli error was a Z error.
void stabilizer_lattice::flip(int x, int y, bool z) {
    if (z) {
        errors_Z[x][y] ^= 1;
        for (const pair<int,int>& coord : offsets_Z) {
            int shift_x = (coord.first + x + w) % w;
            int shift_y = (coord.second + y + h) % h;
            syndromes[shift_x][shift_y] ^= 1;
            clustering_helper.flip(shift_x, shift_y);
        }
    } else {
        errors_not_Z[x][y] ^= 1;
        for (const pair<int,int>& coord : offsets_not_Z) {
            int shift_x = (coord.first + x + w) % w;
            int shift_y = (coord.second + y + h) % h;
            syndromes[shift_x][shift_y] ^= 1;
            clustering_helper.flip(shift_x, shift_y);
        }
    }
}

//checks if the error is correctable, by checking parity on each colored hexagon
bool stabilizer_lattice::correctable(vector<pair<int, int> > &cluster, vector<int> &rect) {
    int colors[3] = {0,0,0};
    for (const pair<int,int>& coord : cluster) {
        int x = coord.first;
        int y = coord.second;
        int color = (x + y) % 3;
        colors[color]++;
    }
    
    return (colors[0] % 2 == colors[1] % 2) && (colors[1] % 2 == colors[2] % 2);
}

//apply a Pauli that removes the error syndrome
void stabilizer_lattice::attempt_correction(vector<pair<int,int>>& cluster, vector<int> &rect) {
    if (cluster.size() > 0) {
        no_bias_correction(cluster, rect);
    }
}

//move an error syndrome to the left (to x_dest)
pair<int,int> stabilizer_lattice::string_clean_x(int x_src, int y_src, int x_dest) {
    int current_x = x_src;
    int current_y = y_src;
    int step = 0;
    
    boost::multi_array<int, 2>* correction_1 = &(correction_Z);
    boost::multi_array<int, 2>* correction_2 = &(correction_not_Z);
    if (stabilizer == 'B') {
        correction_1 = &(correction_not_Z);
        correction_2 = &(correction_Z);
    }
    
    while (current_x != x_dest && current_x != (x_dest - 1 + w) % w) {
        int x1 = (current_x - 1 + w) % w;
        if (step % 2 == 0) {
            (*correction_1)[x1][y_src] ^= 1;
            (*correction_2)[x1][y_src] ^= 1;
            current_x = (current_x - 2 + w) % w;
            current_y = (current_y - 1 + h) % h;
        } else {
            (*correction_1)[current_x][y_src] ^= 1;
            (*correction_2)[x1][y_src] ^= 1;
            current_x = x1;
            current_y = (current_y + 1) % h;
        }
        step++;
    }
    
    return pair<int,int>(current_x,current_y);
}

//move an error syndrome up (to y_dest)
void stabilizer_lattice::string_clean_y(int x_src, int y_src, int x_dest, int y_dest) {
    int current_y = y_src;
    int step = 0;
    int x_dest1 = (x_dest - 1 + w) % w;
    if (x_src == x_dest1) {
        step = 1;
    }
    boost::multi_array<int, 2>* correction_1 = &(correction_Z);
    boost::multi_array<int, 2>* correction_2 = &(correction_not_Z);
    if (stabilizer == 'B') {
        correction_1 = &(correction_not_Z);
        correction_2 = &(correction_Z);
    }
    
    if (y_dest == (y_src + 2) % h) {
        (*correction_2)[x_dest1][(y_src + 1) % h] ^= 1;
        (*correction_1)[x_dest][y_dest] ^= 1;
    } else {
        while (current_y != y_dest) {
            if (step % 2 == 0) {
                (*correction_1)[x_dest][current_y] ^= 1;
                (*correction_2)[x_dest1][(current_y - 1 + h) % h] ^= 1;
                current_y = (current_y - 2 + h) % h;
            } else {
                (*correction_2)[x_dest1][current_y] ^= 1;
                (*correction_1)[x_dest][current_y] ^= 1;
                current_y = (current_y - 1 + h) % h;
            }
            step++;
        }
    }
}

//move an error to point x_dest,y_dest using a stringlike operator
void stabilizer_lattice::string_clean(int x_src, int y_src, int x_dest, int y_dest) {
    pair<int,int> x_cleaned = string_clean_x(x_src, y_src, x_dest);
    string_clean_y(x_cleaned.first, x_cleaned.second, x_dest, y_dest);
}

//correct an error with stringlike operators by moving all errors to the top right of the rect
void stabilizer_lattice::no_bias_correction(std::vector<std::pair<int, int> > &cluster, std::vector<int> &rect) {
    pair<int,int> dests[3];
    int tl_x = rect[0];
    int tl_y = rect[1];
    int tl_y1 = (tl_y + 1 + h) % h;
    int tl_x1 = (tl_x + 1 + w) % w;
    pair<int,int> dests_unordered[3] = {pair<int,int>(tl_x,tl_y),pair<int,int>(tl_x,tl_y1),pair<int,int>(tl_x1,tl_y1)};
    for (const pair<int,int>& coord : dests_unordered) {
        int color = (coord.first + coord.second) % 3;
        dests[color] = coord;
    }
    
    int color_parities[3] = {0,0,0};
    for (const pair<int, int>& syndrome: cluster) {
        int color = (syndrome.first + syndrome.second) % 3;
        string_clean(syndrome.first, syndrome.second, dests[color].first, dests[color].second);
        color_parities[color]++;
    }
    
    if (color_parities[0] % 2 == 1) {
        if (stabilizer == 'A') {
            correction_not_Z[tl_x][tl_y1] ^= 1;
        } else {
            correction_Z[tl_x][tl_y1] ^= 1;
        }
    }
}

/*
 * Try to correct the error. First get clusters. Check if each cluster is correctable, then try to
 * correct each cluster. Repeat with larger cluster sizes until all errors are corrected, or if
 * give_up is true, clusters become larger than min(w,h)/2. In that case return false, else true.
 */
bool stabilizer_lattice::make_correction() {
    fill(correction_Z.data(), correction_Z.data()+correction_Z.num_elements(), 0);
    fill(correction_not_Z.data(), correction_not_Z.data()+correction_not_Z.num_elements(), 0);
    
    int p_M = floor(log2((float)w / 2));
    if (!give_up) {
        p_M++;
    }
    vector<vector<pair<int,int>>> clusters;
    vector<vector<int>> rects;
    int syndromes_remaining = clustering_helper.get_n_syndromes();
    clustering_helper.unmark_all();

    for (int l = 0; l <= p_M; l++) {
        clusters.clear();
        rects.clear();
        clustering_helper.get_clusters(clusters, rects, l);
        
        for (int i = 0; i < clusters.size(); i++) {
            if (give_up && (rects[i][2] >= w/2 || rects[i][3] >= h/2)) {
                return false;
            }
            
            if (correctable(clusters[i], rects[i])) {
                for (const pair<int,int>& coord : clusters[i]) {
                    syndromes_remaining--;
                    clustering_helper.mark_clustered(coord.first, coord.second);
                }
                attempt_correction(clusters[i], rects[i]);
            }
        }
        
        if (syndromes_remaining == 0) {
            if (reduce_weight) {
                reduce_stabilizer_weight();
            }
            return true;
        }
    }
    
    return false;
}

// optional: try to make the correction operator have smaller weight by acting with stabilizers.
void stabilizer_lattice::reduce_stabilizer_weight() {
    int n_passes = 3;
    boost::multi_array<int, 2>* correction_1 = &(correction_Z);
    boost::multi_array<int, 2>* correction_2 = &(correction_not_Z);
    if (stabilizer == 'B') {
        correction_1 = &(correction_not_Z);
        correction_2 = &(correction_Z);
    }
    
    for (int pass = 0; pass < n_passes; pass++) {
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                int count = 0;
                int y1 = (y+1)%h;
                int xp1 = (x+1)%w;
                int xm1 = (x-1+w)%w;
                count += (*correction_1)[x][y];
                count += (*correction_1)[x][y1];
                count += (*correction_1)[xp1][y1];
                count += (*correction_2)[x][y];
                count += (*correction_2)[x][y1];
                count += (*correction_2)[xm1][y];
                
                if (count > 3) {
                    (*correction_1)[x][y] ^= 1;
                    (*correction_1)[x][y1] ^= 1;
                    (*correction_1)[xp1][y1] ^= 1;
                    (*correction_2)[x][y] ^= 1;
                    (*correction_2)[x][y1] ^= 1;
                    (*correction_2)[xm1][y] ^= 1;
                }
            }
        }
    }
}

//for debugging, check if any stabilizers are violated.
bool stabilizer_lattice::check_stabilizers() {
    boost::multi_array<int, 2>* correction_1 = &(correction_Z);
    boost::multi_array<int, 2>* correction_2 = &(correction_not_Z);
    boost::multi_array<int, 2>* errors_1 = &(errors_Z);
    boost::multi_array<int, 2>* errors_2 = &(errors_not_Z);
    if (stabilizer == 'B') {
        correction_1 = &(correction_not_Z);
        correction_2 = &(correction_Z);
        errors_1 = &(errors_not_Z);
        errors_2 = &(errors_Z);
    }
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            int yp1 = (y + 1) % h;
            int xp1 = (x + 1) % w;
            int xm1 = (x - 1 + w) % w;
            
            int parity = 0;
            parity ^= ((*errors_1)[x][y] ^ (*correction_1)[x][y]);
            parity ^= ((*errors_1)[x][yp1] ^ (*correction_1)[x][yp1]);
            parity ^= ((*errors_1)[xp1][yp1] ^ (*correction_1)[xp1][yp1]);
            parity ^= ((*errors_2)[x][y] ^ (*correction_2)[x][y]);
            parity ^= ((*errors_2)[xm1][y] ^ (*correction_2)[xm1][y]);
            parity ^= ((*errors_2)[x][yp1] ^ (*correction_2)[x][yp1]);
            
            if (parity != 0) {
                cout << x << "," << y << " stabilizer failed" << endl;
                return false;
            }
        }
    }
    
    return true;
}

//check if any logical errors are violated. If so correction failed.
bool stabilizer_lattice::check_logicals() {
    boost::multi_array<int, 2>* correction_1 = &(correction_Z);
    boost::multi_array<int, 2>* correction_2 = &(correction_not_Z);
    boost::multi_array<int, 2>* errors_1 = &(errors_Z);
    boost::multi_array<int, 2>* errors_2 = &(errors_not_Z);
    if (stabilizer == 'B') {
        correction_1 = &(correction_not_Z);
        correction_2 = &(correction_Z);
        errors_1 = &(errors_not_Z);
        errors_2 = &(errors_Z);
    }
    
    int red_v = 0;
    int blue_v = 0;
    for (int y = 0; y < h / 3; y++) {
        blue_v ^= ((*correction_1)[1][3*y] ^ (*errors_1)[1][3*y]);
        blue_v ^= ((*correction_1)[1][3*y+2] ^ (*errors_1)[1][3*y+2]);
        blue_v ^= ((*correction_2)[0][3*y] ^ (*errors_2)[0][3*y]);
        blue_v ^= ((*correction_2)[0][3*y+1] ^ (*errors_2)[0][3*y+1]);
        
        red_v ^= ((*correction_1)[1][3*y] ^ (*errors_1)[1][3*y]);
        red_v ^= ((*correction_1)[1][3*y+1] ^ (*errors_1)[1][3*y+1]);
        red_v ^= ((*correction_2)[0][3*y+1] ^ (*errors_2)[0][3*y+1]);
        red_v ^= ((*correction_2)[0][3*y+2] ^ (*errors_2)[0][3*y+2]);
    }
    
    int red_h = 0;
    int blue_h = 0;
    for (int x = 0; x < w / 3; x++) {
        blue_h ^= ((*correction_1)[3*x][0] ^ (*errors_1)[3*x][0]);
        blue_h ^= ((*correction_1)[3*x+2][0] ^ (*errors_1)[3*x+2][0]);
        blue_h ^= ((*correction_2)[3*x][0] ^ (*errors_2)[3*x][0]);
        blue_h ^= ((*correction_2)[3*x+1][0] ^ (*errors_2)[3*x+1][0]);
        
        red_h ^= ((*correction_1)[3*x][0] ^ (*errors_1)[3*x][0]);
        red_h ^= ((*correction_1)[3*x+1][0] ^ (*errors_1)[3*x+1][0]);
        red_h ^= ((*correction_2)[3*x+1][0] ^ (*errors_2)[3*x+1][0]);
        red_h ^= ((*correction_2)[3*x+2][0] ^ (*errors_2)[3*x+2][0]);
    }
    
    return (red_h == 0 && blue_h == 0 && red_v == 0 && blue_v == 0);
}

//check if the correction did not violate stabilizers or logicals
bool stabilizer_lattice::check_correction() {
    return check_logicals();
}

//reset the qubits on the lattice so a new simulation can be performed
void stabilizer_lattice::clear() {
    fill(errors_Z.data(), errors_Z.data()+errors_Z.num_elements(), 0);
    fill(errors_not_Z.data(), errors_not_Z.data()+errors_not_Z.num_elements(), 0);
    fill(correction_Z.data(), correction_Z.data()+correction_Z.num_elements(), 0);
    fill(correction_not_Z.data(), correction_not_Z.data()+correction_not_Z.num_elements(), 0);
    fill(syndromes.data(), syndromes.data()+syndromes.num_elements(), 0);
    clustering_helper.clear();
}
