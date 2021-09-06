//
//  main.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

#include <iostream>
#include <chrono>
#include "threshold_sim.hpp"
#include "thermal_bias_sim.hpp"
#include "ca_sim.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    string filename = "";
    int samples = atoi(argv[1]);
    if (argc == 3) {
        filename = argv[2];
    }

    //system sizes
    vector<pair<int,int>> dims {pair<int,int>(12,15),pair<int,int>(24,27),pair<int,int>(48,51)};
    //time intervals to check for logical errors (for each size)
    vector<int> strides {1,1,1};
    //temperatures
    vector<float> ps {0.0001};
    //bias
    vector<float> normalized_bias {0.0,0.0,1.0};
    //number of steps between using the CA rule
    int ca_freq = 1;

    auto t1 = chrono::high_resolution_clock::now();
    auto timenow =
          chrono::system_clock::to_time_t(chrono::system_clock::now());
    cout << ctime(&timenow) << endl;

    ca_sim::lifetime_sim(samples, dims, strides, ps, normalized_bias, ca_freq, filename);

    auto t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> dur = std::chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    cout << dur.count() << endl;

    return 0;
}
