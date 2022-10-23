//
//  ca_sim.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/27/21.
//

#include "ca_sim.hpp"
#include <iostream>
#include <fstream>

using namespace std;

//w,h = dimensions, bias = prob of I,X,Y,Z, ca_freq = how many steps between using CA rule
ca_sim::ca_sim(int w, int h, vector<float>& bias, int ca_freq) :
w(w),
h(h),
ca_freq(ca_freq),
decoder(w,h,bias),
syndromes_A(boost::extents[w][h]),
syndromes_B(boost::extents[w][h]),
bias_cumulative(4),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
    syndrome_offsets_1 = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    syndrome_offsets_2 = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(1,0)};
    
    bias_cumulative[0] = bias[0];
    for (int i = 1; i < 4; i++) {
        bias_cumulative[i] = bias[i] + bias_cumulative[i-1];
    }
}

//act on qubit at site x,y on sublattice 1/2, with Pauli p
void ca_sim::flip(int x, int y, int sublattice, Pauli p) {
    vector<pair<int,int>>* syndrome_offsets = &syndrome_offsets_1;
    bool flip_A = false;
    bool flip_B = false;
    
    if (sublattice == 1) {
        flip_A = (p == Z || p == Y);
        flip_B = (p == X || p == Y);
    } else {
        flip_A = (p == Y || p == X);
        flip_B = (p == Z || p == X);
        syndrome_offsets = &syndrome_offsets_2;
    }
    
    for (const pair<int,int>& coord : *syndrome_offsets) {
        int shift_x = (x + coord.first + w) % w;
        int shift_y = (y + coord.second + h) % h;
        
        syndromes_A[shift_x][shift_y] ^= flip_A;
        syndromes_B[shift_x][shift_y] ^= flip_B;
    }
    
    decoder.flip(x, y, sublattice, p);
}

//check if the error was correctable by RG decoder and correction+error was trivial
bool ca_sim::check() {
    return !decoder.make_correction() || decoder.check_correction();
}

//adds errors on entire system using single-site probabilities given by bias
void ca_sim::add_errors() {
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            for (int sublattice = 1; sublattice <= 2; sublattice++) {
                float p_flip = dis(engine);
                if (bias_cumulative[0] < p_flip && p_flip <= bias_cumulative[1]) {
                    flip(x, y, sublattice, X);
                } else if (bias_cumulative[1] < p_flip && p_flip <= bias_cumulative[2]) {
                    flip(x, y, sublattice, Y);
                } else if (bias_cumulative[2] < p_flip) {
                    flip(x, y, sublattice, Z);
                }
            }
        }
    }
}

//checks each qubit to see if number of syndromes could be reduced by applying Z on
//the qubit. Threshold = 2 or 3, how many syndromes a Z operator must remove for it to
//be applied.
void ca_sim::ca_correction(int threshold, int sublattice) {
    vector<pair<int,int>>* syndrome_offsets = &syndrome_offsets_1;
    boost::multi_array<int,2>* syndromes = &syndromes_A;
    
    if (sublattice == 2) {
        syndrome_offsets = &syndrome_offsets_2;
        syndromes = &syndromes_B;
    }
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            int syndromes_count = 0;
            
            for (const pair<int, int>& coord : *syndrome_offsets) {
                int shift_x = (x + coord.first + w) % w;
                int shift_y = (y + coord.second + h) % h;
                syndromes_count += (*syndromes)[shift_x][shift_y];
            }
            
            if (syndromes_count >= threshold) {
                flip(x, y, sublattice, Z);
            }
        }
    }
}

//repeat n times, applying CA correction every ca_freq steps
void ca_sim::steps(int n_steps) {
    for (int i = 0; i < n_steps; i++) {
        total_steps++;
        add_errors();
        
        if (total_steps % ca_freq == 0) {
            ca_correction(2, 1);
            ca_correction(2, 2);
        }
    }
}

void ca_sim::debug() {
    decoder.print_error();
    cout << "===========" << endl;
    decoder.print_syndrome();
//    cout << "@@@@@@@@@@@@@@@@@@" << endl;
//    int count = 0;
//    for (int y = 0; y < h; y++) {
//        for (int x = 0; x < w; x++) {
//            cout << syndromes_A[x][y] << " ";
//            count += syndromes_A[x][y];
//        }
//        cout << endl;
//    }
//    cout << endl << count << endl;
}

/*
 Runs a simulation to find lifetime with cellular automaton rule. Each step of the simulation act
 with error channel on all sites, then apply the CA rule every freq steps. Output # of steps until
 logical failure.
 **avg: number of samples
 **dims: array of w,h pairs
 **strides: array of int t, where t is the interval at which to check for logical errors. One per size
 **ps: array of single site error probabilities p = px+py+pz
 **normalized_bias: array of single site error probabilities {px,py,pz} normalized to 1
 **freq: how many steps between each action of the CA rule.
 **filename: file to output to. use "" to not output to file
 */
void ca_sim::lifetime_sim(int avg, vector<pair<int, int> > &dims, vector<int> &strides, vector<float> &ps, vector<float> &normalized_bias, int freq, string filename) {
    cout << "p0,px,py,pz,freq,w,h,t" << endl;
    ofstream myfile;
    if (filename != "") {
        myfile.open(filename);
        myfile << "p0,px,py,pz,freq,w,h,t" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        int stride = strides[i];
        
        for (float p : ps) {
            for (int i = 0; i < avg; i++) {
                vector<float> bias(4);
                bias[0] = 1 - p;
                for (int i = 0; i < 3; i++) {
                    bias[i+1] = normalized_bias[i] * p;
                }
                
                string error_probs = "";
                for (int i = 0; i < 4; i++) {
                    error_probs += to_string(bias[i]) + ",";
                }
                
                cout << error_probs << freq << "," << w << "," << h << ",";
                if (filename != "") {
                    myfile << error_probs << freq << "," << w << "," << h << ",";
                }
                
                ca_sim ca(w,h,bias,freq);
                int t = 0;
                while(ca.check()) {
                    ca.steps(stride);
                    t += stride;
                }
                
                cout << t << endl;
                if (filename != "") {
                    myfile << t << endl;
                }
            }
        }
    }
    
    if (filename != "") {
        myfile.close();
    }
}

void ca_sim::distribution_sim(int avg, int t, std::vector<std::pair<int, int> > &dims, std::vector<float> &ps, std::vector<float> &normalized_bias, int freq, std::string filename) {
    cout << "p0,px,py,pz,freq,w,h,count" << endl;
    ofstream myfile;
    if (filename != "") {
        myfile.open(filename);
        myfile << "p0,px,py,pz,freq,w,h,count" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        
        for (float p : ps) {
            vector<float> bias(4);
            bias[0] = 1 - p;
            for (int i = 0; i < 3; i++) {
                bias[i+1] = normalized_bias[i] * p;
            }
            
            string error_probs = "";
            for (int i = 0; i < 4; i++) {
                error_probs += to_string(bias[i]) + ",";
            }
            
            ca_sim ca(w,h,bias,freq);
            
            for (int i = 0; i < avg; i++) {
                cout << error_probs << freq << "," << w << "," << h << ",";
                if (filename != "") {
                    myfile << error_probs << freq << "," << w << "," << h << ",";
                }
                
                ca.steps(t);
                int syndromes_count = 0;
                for (int x = 0; x < w; x++) {
                    for (int y = 0; y < h; y++) {
                        syndromes_count += (ca.syndromes_A[x][y] + ca.syndromes_B[x][y]);
                    }
                }
                
                cout << syndromes_count << endl;
                if (filename != "") {
                    myfile << syndromes_count << endl;
                }
            }
        }
    }
    
    if (filename != "") {
        myfile.close();
    }
}
