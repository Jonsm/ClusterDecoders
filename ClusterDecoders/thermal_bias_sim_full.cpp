//
//  thermal_bias_sim_full.cpp
//  ClusterDecoders
//
//  Created by Jon on 12/1/21.
//

#include <iostream>
#include <fstream>
#include <limits>
#include <signal.h>
#include "thermal_bias_sim_full.hpp"
#include "exact_decoder.hpp"
#include "cluster_decoder.hpp"

using namespace std;

template <typename T>
bool thermal_bias_sim_full<T>::interrupted = false;

/*
 **w,h = dims
 **bias = p_id,px,py,pz
 **no_correct = whether not to apply CA rule (if no_correct is true, just let errors accumulate)
 **no_correct_bias = experimental, try to correct non Z errors with other CA rule
 */
template <typename T>
thermal_bias_sim_full<T>::thermal_bias_sim_full(int w, int h, std::vector<float>& bias, bool no_correct, bool no_correct_bias) :
w(w),
h(h),
decoder(w,h,bias, false),
s1(w,h,1,bias, no_correct, no_correct_bias),
s2(w,h,2,bias, no_correct, no_correct_bias),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
    s1_offsets = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    s2_offsets = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(1,0)};
    R = s1.R_total() + s2.R_total();
    
    for (int i = 0; i < 3; i++) {
        syndromes_list[i].resize(boost::extents[w][h]);
    }
}

//doesnt work for exact decoder
template <typename T>
int thermal_bias_sim_full<T>::num_errors() {
    return decoder.num_errors();
}

//check if correction is trivial
template<>
bool thermal_bias_sim_full<exact_decoder>::check() {
    return decoder.make_correction();
}

template<>
bool thermal_bias_sim_full<cluster_decoder>::check() {
    return decoder.make_correction() && decoder.check_correction();
}

//apply errors for time >= t_evol
template <typename T>
double thermal_bias_sim_full<T>::mc_time(double t_evol) {
    double t_diff = 0;
    
    while (t_diff < t_evol) {
        if (interrupted) {
            return -1;
        }
        t_diff += wait_rand_time();
        make_flip();
    }
    t += t_diff;
    
    return t;
}

//wait time between flips
template <typename T>
double thermal_bias_sim_full<T>::wait_rand_time() {
    float r = 1.0 - dis(engine);
    while (r == 0) {
        r = 1.0 - dis(engine);
    }
    
    double delta_T = log(r) * -1 / R;
    return delta_T;
}

//flip spin for one sublattice
template <typename T>
void thermal_bias_sim_full<T>::flip_helper(std::pair<int, int> &coord, char (&stabilizers)[2], int sublattice) {
    boost::multi_array<int, 2>* syndromes;
    vector<pair<int, int>>* syndrome_offsets = &s1_offsets;
    if(sublattice == 2) {
        syndrome_offsets = &s2_offsets;
    }
    
    for (int i = 0; i < 2; i++) {
        char stabilizer = stabilizers[i];
        syndromes = &(syndromes_list[stabilizer - 'A']);
        
        for (pair<int,int>& offset : *syndrome_offsets) {
            int shift_x = (offset.first + coord.first + w) % w;
            int shift_y = (offset.second + coord.second + h) % h;
            int syn_old = (*syndromes)[shift_x][shift_y];
            (*syndromes)[shift_x][shift_y] ^= 1;
            
            for (pair<int,int>& offset2 : s1_offsets) {
                int shift_x_2 = (shift_x - offset2.first + w) % w;
                int shift_y_2 = (shift_y - offset2.second + h) % h;
                s1.update_probs(shift_x_2, shift_y_2, stabilizer, syn_old);
            }
            for (pair<int,int>& offset2 : s2_offsets) {
                int shift_x_2 = (shift_x - offset2.first + w) % w;
                int shift_y_2 = (shift_y - offset2.second + h) % h;
                s2.update_probs(shift_x_2, shift_y_2, stabilizer, syn_old);
            }
        }
    }
}

//flip spin
template <typename T>
void thermal_bias_sim_full<T>::make_flip() {
    thermal_sublattice* sublattice = &s1;
    int sublattice_i = 1;
    float sublattice_p = dis(engine) * R;
    if (sublattice_p > s1.R_total()) {
        sublattice = &s2;
        sublattice_i = 2;
    }
    
    Pauli pauli = sublattice->flip_type();
    pair<int,int> coord = sublattice->flip_coord(pauli);
    if (pauli == Y) {
        y_count++;
    }
    
    char stabilizers[2];
    stabilizers[0]='A';
    stabilizers[1]='B';
    if (sublattice_i == 1 && pauli != Y) {
        if (pauli == X) {
            stabilizers[0]='B';
        }
        stabilizers[1]='C';
    } else if (sublattice_i == 2 && pauli != X) {
        if (pauli == Z) {
            stabilizers[0]='B';
        }
        stabilizers[1]='C';
    }
    
    flip_helper(coord, stabilizers, sublattice_i);
    R = s1.R_total() + s2.R_total();
    decoder.flip(coord.first, coord.second, sublattice_i, pauli);
}

//debug
template <typename T>
void thermal_bias_sim_full<T>::debug() {
    cout << endl;
    decoder.print_error();
    cout << "=============" << endl;
    decoder.print_syndrome();
}

//error handler in case interrupted by cluster
template <typename T>
void thermal_bias_sim_full<T>::signal_callback_handler(int signum) {
    interrupted = true;
}

/*
 **average = number of runrs
 **strides = times (for each size) to wait between error checks
 **dims = system sizes
 **ps = total error rates
 **normalized_biases = normalized error rate px,py,pz. Will run entire simulation for
   each combination of normalized_bias, ps
 **filename = output filename
 **no_correct = as in thermal_bias_sim
 **no_correct_bias = as in thermal_bias_sim
 */
template <typename T>
void thermal_bias_sim_full<T>::lifetime_sim(int average, std::vector<double> &strides, std::vector<std::pair<int, int> > &dims, std::vector<float> &ps, std::vector<std::vector<float>> &normalized_biases, std::string filename, bool no_correct, bool no_correct_bias) {
    bool output_to_file = (filename != "");
    
    signal(SIGINT, signal_callback_handler);
    signal(SIGTERM, signal_callback_handler);
    
    cout << "p0,px,py,pz,w,h,t,y_count" << endl;
    ofstream myfile;
    if (output_to_file) {
        myfile.open(filename);
        myfile << "p0,px,py,pz,w,h,t,y_count" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        double stride = strides[i];
        
        for (vector<float>& normalized_bias : normalized_biases) {
            for (float p : ps) {
                vector<float> bias(4);
                bias[0] = 1 - p;
                for (int i = 0; i < 3; i++) {
                    bias[i+1] = normalized_bias[i] * p;
                }
                
                stringstream error_probs_stream;
                for (int i = 0; i < 4; i++) {
                    error_probs_stream << bias[i] << ",";
                }
                string error_probs = error_probs_stream.str();
                
                for (int i = 0; i < average; i++) {
                    cout << error_probs << w << "," << h << ",";
                    if (filename != "") {
                        myfile << error_probs << w << "," << h << ",";
                    }
                    
                    thermal_bias_sim_full<T> sim(w,h,bias, no_correct, no_correct_bias);
                    double t = 0;
                    bool check = true;
                    do {
                        t = sim.mc_time(stride);
                        if (interrupted) {
                            cout << std::numeric_limits<float>::max() << "," << sim.y_count << endl;
                            if (filename != "") {
                                myfile << std::numeric_limits<float>::max() << "," << sim.y_count << endl;
                            }
                            return;
                        }
                        check = sim.check();
                    } while (check);
                    
                    cout << t << "," << sim.y_count << endl;
                    if (filename != "") {
                        myfile << t << "," << sim.y_count << endl;
                    }
                }
            }
        }
    }
    
    if (output_to_file) {
        myfile.close();
    }
}

/*
 **average = number of runs
 **strides = times (for each size) to wait between error checks
 **dims = system sizes
 **ps = total error rates
 **bias_start = smallest value of py/pz
 **bias_stop = largest value of py/pz
 **bias_steps = number of steps
 **filename = output filename
 **no_correct = as in thermal_bias_sim
 **no_correct_bias = as in thermal_bias_sim
 
 * Only real difference here is how it outputs results.
 */
template <>
void thermal_bias_sim_full<cluster_decoder>::compare_bias_sim(int average, std::vector<double> &strides, std::vector<std::pair<int, int>> &dims, std::vector<float> &ps, float bias_start, float bias_stop, int bias_steps, std::string filename, bool no_correct_bias) {
    bool output_to_file = (filename != "");
    
    cout << "CA_on,p,py,w,h,t" << endl;
    ofstream myfile;
    if (output_to_file) {
        myfile.open(filename);
        myfile << "CA_on,p,py,w,h,t" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        double stride = strides[i];
        for (float p : ps) {
            float py_inc = (bias_stop - bias_start) / (bias_steps - 1);
            for (int j = 0; j < bias_steps; j++) {
                float py = p*bias_start + p*py_inc*j;
                vector<float> bias(4);
                bias[0] = 1 - p;
                bias[3] = p - py;
                bias[2] = py;
                
                for (int k = 0; k < average; k++) {
                    for (int no_correct = 0; no_correct <= 1; no_correct++) {
                        cout << !no_correct << "," << p << "," << py << "," << w << "," << h << ",";
                        if (output_to_file) {
                            myfile << !no_correct << "," << p << "," << py << "," << w << "," << h << ",";
                        }
                        
                        thermal_bias_sim_full<cluster_decoder> sim(w,h,bias, no_correct, no_correct_bias);
                        double t = 0;
                        bool check = true;
                        do {
                            t = sim.mc_time(stride);
                            check = sim.check();
                        } while (check);
                        
                        cout << t << endl;
                        if (filename != "") {
                            myfile << t << endl;
                        }
                    }
                }
            }
        }
    }
    
    if (output_to_file) {
        myfile.close();
    }
}

template class thermal_bias_sim_full<exact_decoder>;
template class thermal_bias_sim_full<cluster_decoder>;
