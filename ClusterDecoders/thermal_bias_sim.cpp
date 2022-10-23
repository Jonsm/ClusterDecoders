//
//  thermal_bias_sim.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/22/21.
//

//DEPRECATED

#include "thermal_bias_sim.hpp"
#include "cluster_decoder.hpp"
#include <iostream>
#include <fstream>

using namespace std;

//create array of transition rates for qubits, depending on which neighbors are flipped
vector<float> get_unique_gammas(float T) {
    std::vector<float> unique_gammas;
    double beta = 1.0 / T;
    for (double omega = -3; omega <= 3; omega += 2) {
        double gamma = omega / (1.0 - exp(-beta * omega));
        unique_gammas.push_back(gamma);
    }
    return unique_gammas;
}

//w,h = dims, T = temperature, stabilizer = which stabilizer to simulate. shares decoder with
//a thermal_bias_sim for the other stabilizer
thermal_bias_sim::thermal_bias_sim(int w, int h, float T, char stabilizer, cluster_decoder* decoder) :
w(w),
h(h),
T(T),
stabilizer(stabilizer),
decoder(decoder),
syndromes(boost::extents[w][h]),
gamma_grid(boost::extents[w][h]),
rd_array(get_unique_gammas(T), w*h),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
    if (stabilizer == 'A') {
        syndrome_offsets  = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    } else {
        syndrome_offsets  = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(1,0)};
    }
    
    R = rd_array.total_prob();
}

//act with the error channel for a time no smaller than t_evol. Returns the time of
//the first flip after t_evol
double thermal_bias_sim::mc_time(double t_evol) {
    double t_diff = 0;
    
    while (t_diff < t_evol) {
        t_diff += wait_rand_time();
        make_flip();
    }
    t += t_diff;
    
    return t;
}

//wait a time given by exponential distribution
double thermal_bias_sim::wait_rand_time() {
    float r = 1.0 - dis(engine);
    while (r == 0) {
        r = 1.0 - dis(engine);
    }
    
    double delta_T = log(r) * -1 / R;
    return delta_T;
}

//flip a random qubit, where qubits are more likely to be flipped if energy is reduced
void thermal_bias_sim::make_flip() {
    pair<int,int> coord = rand_site();
    int x = coord.first;
    int y = coord.second;
    
    for (const pair<int, int>& coord2 : syndrome_offsets) {
        int syndrome_before = syndromes[(x + coord2.first + w) % w][(y + coord2.second + h) % h];
        int x2 = (x + coord2.first + w) % w;
        int y2 = (y + coord2.second + h) % h;
        syndromes[x2][y2] ^= 1;
        
        for (const pair<int, int>& coord3 : syndrome_offsets) {
            int x3 = (x2 - coord3.first + w) % w;
            int y3 = (y2 - coord3.second + h) % h;
            int gg_before = gamma_grid[x3][y3];
            pair<int, int> coord4(x3,y3);
            
            int gg_new = gg_before + (1 - 2 * syndrome_before);
            gamma_grid[x3][y3] = gg_new;
            int ind = y3 * w + x3;
            rd_array.assign(ind, gg_new);
        }
    }
    R = rd_array.total_prob();
    
    int sublattice = 1;
    if (stabilizer == 'B') {
        sublattice = 2;
    }
    decoder->flip(x,y,sublattice,Z);
}

//find a random site, where the site is more likely to be picked if flipping it reduces energy
pair<int,int> thermal_bias_sim::rand_site() {
    int ind = rd_array.get_random();
    int x = ind % w;
    int y = ind / w;
    return pair<int, int>(x,y);
}

/*
 Runs a simulation to find lifetime with thermal Z noise. Method from sec. 3 of https://arxiv.org/pdf/1411.6643.pdf. Outputs first time of logical failure.
 **avg: number of samples
 **dims: array of w,h pairs
 **strides: array of float t, where t is the interval at which to check for logical errors. One per size
 **Ts: array of temperatures
 **filename: file to output to. use "" to not output to file
 */
void thermal_bias_sim::lifetime_sim(int average, std::vector<double> &strides, std::vector<std::pair<int, int> > &dims, std::vector<float> &Ts, std::string filename) {
    bool output_to_file = (filename != "");
    
    cout << "T,w,h,t" << endl;
    ofstream myfile;
    if (output_to_file) {
        myfile.open(filename);
        myfile << "T,w,h,t" << endl;
    }
    
    for (int i = 0; i < dims.size(); i++) {
        int w = dims[i].first;
        int h = dims[i].second;
        double stride = strides[i];
        
        for (float T : Ts) {
            for (int i = 0; i < average; i++) {
                cout << T << "," << w << "," << h << ",";
                if (output_to_file) {
                    myfile << T << "," << w << "," << h << ",";
                }
                    
                vector<float> bias {1.0,0.0,0.0,0.0};
                cluster_decoder decoder(w, h, bias);
                thermal_bias_sim sim_A(w,h,T,'A', &decoder);
                thermal_bias_sim sim_B(w,h,T,'B', &decoder);
                
                double t_A = 0, t_B = 0;
                do {
                    t_A = sim_A.mc_time(stride);
                    t_B = sim_B.mc_time(stride);
                } while (!decoder.make_correction() || decoder.check_correction());
                
                double t = min(t_A, t_B);
                cout << t << endl;
                if (output_to_file) {
                    myfile << t << endl;
                }
            }
        }
    }
    
    if (output_to_file) {
        myfile.close();
    }
}
