//
//  thermal_bias_sim_full.hpp
//  ClusterDecoders
//
//  Created by Jon on 12/1/21.
//

/*
 * Helper class for testing the decoders. Runs BKL algorithm. Runs either:
 ** Lifetime sim: Get memory time for a range of system sizes and error channels.
 ** Compare bias sim: similar, but sweeps bias between a start and stop value.
 */

#ifndef thermal_bias_sim_full_hpp
#define thermal_bias_sim_full_hpp

#include "thermal_sublattice.hpp"
#include "boost/multi_array.hpp"

template <typename T>
class thermal_bias_sim_full {
public:
    int y_count=0;
    
    thermal_bias_sim_full(int w, int h, std::vector<float>& bias, bool no_correct=false, bool no_correct_bias=false);
    double mc_time(double t_evol);
    bool check();
    void debug();
    int num_errors();
    
    static bool interrupted;
    static void lifetime_sim(int average, std::vector<double>& strides, std::vector<std::pair<int,int>>& dims, std::vector<float>& ps, std::vector<std::vector<float>>& normalized_biases, std::string filename, bool no_correct=false, bool no_correct_bias=false);
    static void compare_bias_sim(int average, std::vector<double>& strides, std::vector<std::pair<int,int>>& dims, std::vector<float>& ps, float bias_start, float bias_stop, int bias_steps, std::string filename, bool no_correct_bias=false);
    
private:
    double t = 0;
    float R;
    const int w;
    const int h;
    boost::multi_array<int, 2> syndromes_list[3];
    std::vector<std::pair<int,int>> s1_offsets;
    std::vector<std::pair<int,int>> s2_offsets;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    thermal_sublattice s1;
    thermal_sublattice s2;
    T decoder;
    
    double wait_rand_time();
    void flip_helper(std::pair<int,int>& coord, char (&stabilizers)[2], int sublattice);
    void make_flip();
    
    static void signal_callback_handler(int signum);
};

#endif /* thermal_bias_sim_full_hpp */
