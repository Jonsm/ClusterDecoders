//
//  random_ca_sim.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 10/26/21.
//

//DEPRECATED

#include <iostream>
#include <fstream>
#include "math.h"
#include "random_ca_sim.hpp"

using namespace std;

random_ca_sim::random_ca_sim(int w, int h, vector<float>& bias) :
ca_sim(w, h, bias, 1),
fs(4),
p(bias[3])
{
    init_constants();
}

void random_ca_sim::init_constants() {
    beta = -7.0*log(p)/24;
    dt = (2.0*p/3)*(exp(beta * 3) - 1);
    
    vector<float> Gs;
    for (int nj = -3; nj <= 3; nj+=2) {
        float G_nj = dt * nj / (1.0 - exp(beta * nj * -1));
        Gs.push_back(G_nj);
        cout << "G:" << G_nj << endl;
    }
    
    for (int i = 0; i < 4; i++) {
        float G_nj = Gs[i];
        float G_nj_m = Gs[3-i];
        fs[i] = ((1.0 - p)*G_nj + p*G_nj_m - p) / (1.0 - 2.0*p);
        cout << "F:" << fs[i] << endl;
    }
    cout << "beta:" << beta << endl;
    cout << "dt:" << dt << endl;
}

void random_ca_sim::ca_correction(int threshold, int sublattice) {
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
            
            float f = fs[syndromes_count];
            if (dis(engine) < f) {
                flip(x, y, sublattice, Z);
            }
        }
    }
}

void random_ca_sim::distribution_sim(int avg, int t, std::vector<std::pair<int, int> > &dims, std::vector<float> &ps, std::vector<float> &normalized_bias, std::string filename) {
    cout << "p0,px,py,pz,freq,w,h,count" << endl;
    ofstream myfile;
    if (filename != "") {
        myfile.open(filename);
        myfile << "p0,px,py,pz,w,h,count" << endl;
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
            
            random_ca_sim ca(w,h,bias);
            
            for (int i = 0; i < avg; i++) {
                cout << error_probs << w << "," << h << ",";
                if (filename != "") {
                    myfile << error_probs << w << "," << h << ",";
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
            
            ca.debug();
        }
    }
    
    if (filename != "") {
        myfile.close();
    }
}
