//
//  thermal_bias_sim.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/22/21.
//

// thermal simulation with the spin boson noise model at 100% bias.

#ifndef thermal_bias_sim_hpp
#define thermal_bias_sim_hpp

#include <vector>
#include "random_discrete_array.hpp"
#include "boost/multi_array.hpp"

class cluster_decoder;

class thermal_bias_sim {
public:
    thermal_bias_sim(int w, int h, float T, char stabilizer, cluster_decoder* decoder);
    double mc_time(double t_evol);
    
    //run a thermal simulation and print to file
    static void lifetime_sim(int average, std::vector<double>& strides, std::vector<std::pair<int,int>>& dims, std::vector<float>& Ts, std::string filename);
private:
    double t = 0;
    float R;
    const float T;
    const int w;
    const int h;
    const char stabilizer;
    boost::multi_array<int, 2> syndromes;
    boost::multi_array<int, 2> gamma_grid;
    std::vector<std::pair<int,int>> syndrome_offsets;
    std::vector<float> unique_gammas;
    random_discrete_array rd_array;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    cluster_decoder* decoder;
    
    double wait_rand_time();
    void make_flip();
    std::pair<int,int> rand_site();
};

#endif /* thermal_bias_sim_hpp */
