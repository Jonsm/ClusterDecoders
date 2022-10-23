//
//  thermal_finite_sim.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 11/1/21.
//

//DEPRECATED

#include <iostream>
#include <fstream>
#include <sstream>
#include "math.h"
#include "thermal_finite_sim.hpp"
#include "cluster_decoder.hpp"

using namespace std;

float thermal_finite_sim::get_beta(std::vector<float> &bias) {
    float p_z = bias[3];
    return -7.0*log(p_z)/24;
}

thermal_finite_sim::thermal_finite_sim(int w, int h, std::vector<float>& bias) :
w(w),
h(h),
decoder(w, h, bias, true),
stabs_A(w, h, 1.0/get_beta(bias), 'A'),
stabs_B(w, h, 1.0/get_beta(bias), 'B'),
beta(get_beta(bias)),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
    dt = (2.0*bias[3]/3)*(exp(beta * 3) - 1);
//    cout << "beta" << beta << ",dt" << dt << endl;
    
    float omega_X = bias[1] / dt;
    float omega_Y = bias[2] / dt;
    R_X_total = w*h*2*omega_X;
    R_Y_total = w*h*2*omega_Y;
}

double thermal_finite_sim::mc_time(double t_evol) {
    double t_diff = 0;
    
    float R_avg = 0;
    int n = 0;
    
    while (t_diff < t_evol) {
        t_diff += wait_rand_time();
        make_flip();
        n++;
        R_avg += (R_Y_total / R);
    }
    t += t_diff;
//    cout << "#######" << t << endl;
//    cout << n << "------" << R_avg << "----------" << R_Y_total << endl;
    
    return t;
}

long thermal_finite_sim::mc_timesteps(long steps) {
    double t = steps * dt;
    double t_actual = mc_time(t);
    return (long)(t_actual / dt);
}

bool thermal_finite_sim::check() {
    return decoder.make_correction() && decoder.check_correction();
}

void thermal_finite_sim::debug() {
//    decoder.print_error();
    cout << endl;
//    decoder.print_syndrome();
    cout << "===================" << endl;
    stabs_A.debug();
    stabs_B.debug();
}

void thermal_finite_sim::make_flip() {
    float p_flip_type = dis(engine) * R;
    if (p_flip_type <= R - R_X_total - R_Y_total) {
        if (p_flip_type < stabs_A.get_R()) {
            pair<int,int> coord = stabs_A.rand_site_thermal();
            stabs_A.make_flip(coord.first, coord.second, true);
            decoder.flip(coord.first,coord.second, 1, Z);
        } else {
            pair<int,int> coord = stabs_B.rand_site_thermal();
            stabs_B.make_flip(coord.first, coord.second, true);
            decoder.flip(coord.first,coord.second, 2, Z);
        }
    } else {
        int sublattice = rand() % 2 + 1;
        int x = rand() % w;
        int y = rand() % h;
        if (p_flip_type <= R - R_Y_total) {
            if (sublattice == 1) {
                stabs_B.make_flip(x, y, false);
            } else {
                stabs_A.make_flip(x, y, false);
                stabs_B.make_flip(x, y, true);
            }
            decoder.flip(x, y, sublattice, X);
        } else {
//            cout << "YFLIP:" << p_flip_type << "," << R << "," << R_Y_total << endl;
            if (sublattice == 1) {
                stabs_B.make_flip(x, y, false);
                stabs_A.make_flip(x, y, true);
            } else {
                stabs_A.make_flip(x, y, false);
            }
            decoder.flip(x, y, sublattice, Y);
        }
    }
}

double thermal_finite_sim::wait_rand_time() {
    float r = 1.0 - dis(engine);
    while (r == 0) {
        r = 1.0 - dis(engine);
    }
    
    R = stabs_A.get_R() + stabs_B.get_R() + R_X_total + R_Y_total;
//    cout << stabs_A.get_R() + stabs_B.get_R() << "," << R_Y_total << endl;
    double delta_T = log(r) * -1 / R;
    return delta_T;
}

void thermal_finite_sim::lifetime_sim(int average, std::vector<long> &strides, std::vector<std::pair<int, int> > &dims, std::vector<float> &ps, std::vector<std::vector<float>> &normalized_biases, std::string filename) {
    bool output_to_file = (filename != "");
    
    cout << "p0,px,py,pz,w,h,t" << endl;
    ofstream myfile;
    if (output_to_file) {
        myfile.open(filename);
        myfile << "p0,px,py,pz,w,h,t" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        long stride = strides[i];
        
        for (vector<float>& normalized_bias : normalized_biases) {
            for (float p : ps) {
                vector<float> bias(4);
                bias[0] = 1 - p;
                for (int i = 0; i < 3; i++) {
                    bias[i+1] = normalized_bias[i] * p;
                }
                
//                string error_probs = "";
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
                    
                    thermal_finite_sim sim(w,h,bias);
                    double t = 0;
                    double stride_time = sim.dt * stride;
        //            cout << sim.dt << endl;
                    while (sim.check()) {
                        t = sim.mc_time(stride_time);
                    }
    //                sim.debug();
                    
                    cout << t << endl;
                    if (filename != "") {
                        myfile << t << endl;
                    }
                }
            }
        }
    }
    
    if (output_to_file) {
        myfile.close();
    }
}
