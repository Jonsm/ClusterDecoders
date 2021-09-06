//
//  threshold_sim.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/13/21.
//

#include <fstream>
#include <iostream>
#include "threshold_sim.hpp"

using namespace std;

/*
 **w,h = dims
 **bias = array of {pI,px,py,pz}
 **give_up = whether the RG decoder should report failure if clusters are larger than min(w,h)/2
 **reduce_weight = whether decoder should reduce weight of correction by multiplying with stabilizers
   (good for debugging)
 */
threshold_sim::threshold_sim(int w, int h, vector<float> bias, bool give_up, bool reduce_weight) :
w(w),
h(h),
decoder(w, h, bias, give_up, reduce_weight),
dis(0.0,1.0),
bias_cumulative(4)
{
    bias_cumulative[0] = bias[0];
    for (int i = 1; i < 4; i++) {
        bias_cumulative[i] = bias[i] + bias_cumulative[i-1];
    }
    engine.seed(random_device{}());
}

//flips qubits at x,y in both sublattices, with error probabilites given by bias
void threshold_sim::single_qubit_channel(int x, int y) {
    for (int sublattice = 1; sublattice <= 2; sublattice++) {
        float p_flip = dis(engine);
        if (bias_cumulative[0] < p_flip && p_flip <= bias_cumulative[1]) {
            decoder.flip(x, y, sublattice, X);
        } else if (bias_cumulative[1] < p_flip && p_flip <= bias_cumulative[2]) {
            decoder.flip(x, y, sublattice, Y);
        } else if (bias_cumulative[2] < p_flip) {
            decoder.flip(x, y, sublattice, Z);
        }
    }
}

//flip all qubits and then check if decoding is successful
bool threshold_sim::single_run() {
    decoder.clear();
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            single_qubit_channel(x, y);
        }
    }

    bool corrected = decoder.make_correction();
    if (corrected) {
        return decoder.check_correction();
    } else {
        return false;
    }
}

void threshold_sim::debug() {
    decoder.print_error();
    decoder.print_syndrome();
    decoder.print_correction();
    cout << "~~~~~~~~~~" << endl;
    decoder.print_error_after_correction();
    cout << "########################" << endl << endl;
}

/*
 Outputs data to determine the threshold of the decoder.
 **dims = array of pairs w,h
 **normalized_bias = array {px,py,pz} normalized to add up to 1
 **p_start,p_stop,p_steps = range and # of steps to check error rate p, where p=px+py+pz
 **samples = number of samples
 **filename = file to output to. Leave "" for no output
 **give_up = whether decoder reports failure if clusters are larger than min(w,h)/2
 **reduce_weight = whether decoder should reduce weight of correction by multiplying with stabilizers
 (good for debugging)
 */
void threshold_sim::run_sim(std::vector<std::pair<int, int> > &dims, std::vector<float> normalized_bias, float p_start, float p_stop, int p_steps, int samples, std::string filename, bool give_up, bool reduce_weight) {
    ofstream file;
    if (filename != "") {
        file.open(filename);
    }
    
    cout << "w,h,p,samples,success" << endl;
    if (filename != "") {
        file << "w,h,p,samples,success" << endl;
    }
    
    float step = (p_stop - p_start) / p_steps;
    for (float p = p_start; p <= p_stop; p += step) {
        vector<float> ps(4);
        ps[0] = 1 - p;
        for (int i = 0; i < 3; i++) {
            ps[i+1] = normalized_bias[i] * p;
        }
        
        for (int i = 0; i < dims.size(); i++) {
            int w = dims[i].first;
            int h = dims[i].second;
            
            threshold_sim sim(w, h, ps, give_up, reduce_weight);
            int n_success = 0;
            for (int i = 0; i < samples; i++) {
                bool success = sim.single_run();
                n_success += success;
            }
            
            cout << w << "," << h << "," << p << "," << samples << "," << n_success << endl;
            if (filename != "") {
                file << w << "," << h << "," << p << "," << samples << "," << n_success << endl;
            }
        }
    }
    
    if (filename != "") {
        file.close();
    }
}
