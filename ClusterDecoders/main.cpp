//
//  main.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

#include "math.h"
#include <iostream>
#include <chrono>
#include "thermal_bias_sim_full.hpp"
#include "exact_decoder.hpp"
#include "cluster_decoder.hpp"
#include "threshold_sim.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    string filename = "";
    int samples = atoi(argv[1]);
    if (argc >= 3) {
        filename = argv[2];
    }

    //system sizes
    vector<pair<int,int>> dims {pair<int,int>(6,9),pair<int,int>(12,15),pair<int,int>(24,27),pair<int,int>(48,51),pair<int,int>(96,99)};
    //time interval, for each system size, at which to check corrction
    vector<double> strides {10,10,10,10,10};
    
    float bias_start = .01;
    float bias_stop = .1;
    int bias_steps = 10;
    
    vector<float> ps {1e-5};
    
    if (argc == 4) {
        float p_min = 1e-6;
        float p_step = 1e-6;
        
        int ind = atoi(argv[3]);
        float p = p_min + p_step * ind;
        ps[0] = p;
    }

    auto t1 = chrono::high_resolution_clock::now();
    auto timenow =
          chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << ctime(&timenow) << endl;
    
    //run simulation
    thermal_bias_sim_full<cluster_decoder>::compare_bias_sim(samples, strides, dims, ps, bias_start, bias_stop, bias_steps, filename, true);

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> dur = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << dur.count() << endl;

    return 0;
}
