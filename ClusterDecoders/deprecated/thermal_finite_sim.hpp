//
//  thermal_finite_sim.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 11/1/21.
//

//DEPRECATED

#ifndef thermal_finite_sim_hpp
#define thermal_finite_sim_hpp

#include <vector>
#include <random>
#include "cluster_decoder.hpp"
#include "thermal_finite_stabilizers.hpp"

class thermal_finite_sim {
public:
    thermal_finite_sim(int w, int h, std::vector<float>& bias);
    double mc_time(double t_evol);
    long mc_timesteps(long steps);
    bool check();
    void debug();
    
    static void lifetime_sim(int average, std::vector<long>& strides, std::vector<std::pair<int,int>>& dims, std::vector<float>& ps, std::vector<std::vector<float>>& normalized_biases, std::string filename);
    
private:
    double t = 0;
    const int w;
    const int h;
    float beta;
    float dt;
    float R;
    float R_X_total;
    float R_Y_total;
    cluster_decoder decoder;
    thermal_finite_stabilizers stabs_A;
    thermal_finite_stabilizers stabs_B;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    
    double wait_rand_time();
    void make_flip();
    
    static float get_beta(std::vector<float>& bias);
};

#endif /* thermal_finite_sim_hpp */
