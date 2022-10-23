//
//  threshold_sim.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/13/21.
//

// simulation for finding the threshold of the code at arbitrary bias

#ifndef threshold_sim_hpp
#define threshold_sim_hpp

#include <vector>
#include <random>
#include "paulis.hpp"

template <typename T>
class threshold_sim {
public:
    threshold_sim(int w, int h, std::vector<float> bias, bool give_up, bool reduce_weight);
    bool single_run();
    void debug();
    
    //run a simulation to find the threshold. normalized_bias is the probabilities of X,Y,Z Pauli errors
    //normalized so their sum is 1. Total probability of error (i.e. P_X+P_Y+P_Z) runs from p_start to
    //p_stop in p_steps. Outputs to file.
    static void run_sim(std::vector<std::pair<int,int>>& dims, std::vector<float> normalized_bias, float p_start, float p_stop, int p_steps, int samples, std::string filename, bool give_up=true, bool reduce_weight=false);
    static void run_bias_sim(std::vector<std::pair<int, int> > &dims, std::vector<float> biases, float p_start, float p_stop, int p_steps, int samples, std::string filename, bool give_up=true, bool reduce_weight=false);
private:
    int w;
    int h;
    std::vector<float> bias_cumulative;
    T decoder;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    
    
    void single_qubit_channel(int x, int y);
};
#endif /* threshold_sim_hpp */
