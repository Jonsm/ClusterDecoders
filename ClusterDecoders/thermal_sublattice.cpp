//
//  thermal_sublattice.cpp
//  ClusterDecoders
//
//  Created by Jon on 12/1/21.
//

#include <iostream>
#include "math.h"
#include "thermal_sublattice.hpp"

using namespace std;

/*
 **w,h = dims
 **sublattice = 1 or 2
 **bias = p_id,px,py,pz
 **no_correct = whether not to apply CA rule (if no_correct is true, just let errors accumulate)
 **no_correct_bias = experimental, try to correct non Z errors with other CA rule
 */
thermal_sublattice::thermal_sublattice(int w, int h, int sublattice, vector<float>& bias, bool no_correct, bool no_correct_bias) :
w(w),
h(h),
sublattice(sublattice),
x_rdarray(get_unique_gammas(bias, X, no_correct, no_correct_bias), w*h),
y_rdarray(get_unique_gammas(bias, Y, no_correct, no_correct_bias), w*h),
z_rdarray(get_unique_gammas(bias, Z, no_correct, no_correct_bias), w*h),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
}

//get which type of error to flip, based on bias
Pauli thermal_sublattice::flip_type() {
    float R = x_rdarray.total_prob() + y_rdarray.total_prob() + z_rdarray.total_prob();
    float flip_prob = dis(engine) * R;
    if (flip_prob < x_rdarray.total_prob()) {
        return X;
    } else if (flip_prob < R - z_rdarray.total_prob()) {
        return Y;
    } else {
        return Z;
    }
}

//find which coord to flip, based on BKL algorithm
pair<int,int> thermal_sublattice::flip_coord(Pauli pauli) {
    int coord_ind = 0;
    if (pauli == X) {
        coord_ind = x_rdarray.get_random();
    } else if (pauli == Y) {
        coord_ind = y_rdarray.get_random();
    } else {
        coord_ind = z_rdarray.get_random();
    }
    
    return pair<int,int>{coord_ind % w, coord_ind / w};
}

//update probabilities of flipping each site after flip
void thermal_sublattice::update_probs(int x, int y, char stabilizer, int old_st) {
    int ind = y * w + x;
    
    random_discrete_array* arrs[2];
    if (sublattice == 1) {
        if (stabilizer == 'A') {
            arrs[0] = &y_rdarray;
            arrs[1] = &z_rdarray;
        } else if (stabilizer == 'B') {
            arrs[0] = &x_rdarray;
            arrs[1] = &y_rdarray;
        } else {
            arrs[0] = &x_rdarray;
            arrs[1] = &z_rdarray;
        }
    } else {
        if (stabilizer == 'A') {
            arrs[0] = &x_rdarray;
            arrs[1] = &y_rdarray;
        } else if (stabilizer == 'B') {
            arrs[0] = &x_rdarray;
            arrs[1] = &z_rdarray;
        } else {
            arrs[0] = &y_rdarray;
            arrs[1] = &z_rdarray;
        }
    }
    
    for (int i = 0; i < 2; i++) {
        random_discrete_array* rd_array = arrs[i];
        int old_p_index = rd_array->get_prob_index(ind);
        int new_p_index = old_p_index + (1 - 2*old_st);
        rd_array->assign(ind, new_p_index);
    }
}

//get total rate of flipping any site
float thermal_sublattice::R_total() {
    return x_rdarray.total_prob() + z_rdarray.total_prob() + y_rdarray.total_prob();
}

//debug
void thermal_sublattice::debug() {
    cout << "SUBLATTICE::" << sublattice << endl;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int ind = y * w + x;
            cout << y_rdarray.get_prob_index(ind) << " ";
        }
        cout << endl;
    }
}

//get gammas (flip rates), for each site, based on the bias and what kind of CA rule is used
vector<float> thermal_sublattice::get_unique_gammas(std::vector<float> &bias, Pauli pauli, bool no_correct, bool no_correct_bias) {
    float pz = bias[3];
    float beta = (1.0/6)*log((6.0+pz)/pz);
    vector<float> gammas(7);
    
    if (no_correct) {
        if (pauli == Z) {
            fill(gammas.begin(), gammas.end(), bias[3]);
        } else if (pauli == X) {
            fill(gammas.begin(), gammas.end(), bias[1]);
        } else {
            fill(gammas.begin(), gammas.end(), bias[2]);
        }
    } else {
        if (no_correct_bias && pauli != Z) {
            if (pauli == X) {
                fill(gammas.begin(), gammas.end(), bias[1]);
            } else if (pauli == Y) {
                fill(gammas.begin(), gammas.end(), bias[2]);
            }
        } else {        
            for (int i = 0; i < 7; i++) {
                int omega = -6 + 2*i;
                
                float gamma = beta;
                if (omega != 0) {
                    gamma = (float)omega / (1.0-exp(-1.0*beta*omega));
                }
                
                if (pauli == X) {
                    gamma *= (bias[1]/pz);
                } else if (pauli == Y) {
                    gamma *= (bias[2]/pz);
                }
                
                gammas[i] = gamma;
            }
        }
    }
    return gammas;
}
