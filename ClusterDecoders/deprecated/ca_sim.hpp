//
//  ca_sim.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/27/21.
//

#ifndef ca_sim_hpp
#define ca_sim_hpp

#include <vector>
#include <random>
#include "boost/multi_array.hpp"
#include "cluster_decoder.hpp"
#include "paulis.hpp"

class ca_sim {
public:
    boost::multi_array<int, 2> syndromes_A;
    boost::multi_array<int, 2> syndromes_B;
    
    ca_sim(int w, int h, std::vector<float>& bias, int ca_freq);
    void steps(int n_steps);
    bool check();
    void debug();
    
    static void lifetime_sim(int avg, std::vector<std::pair<int,int>>& dims, std::vector<int>& strides, std::vector<float>& ps, std::vector<float>& normalized_bias, int freq, std::string filename);
    static void distribution_sim(int avg, int t, std::vector<std::pair<int,int>>& dims, std::vector<float>& ps, std::vector<float>& normalized_bias, int freq, std::string filename);
protected:
    const int w;
    const int h;
    const int ca_freq;
    int total_steps = 0;
    std::vector<float> bias_cumulative;
    std::vector<std::pair<int, int>> syndrome_offsets_1;
    std::vector<std::pair<int, int>> syndrome_offsets_2;
    cluster_decoder decoder;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    
    void flip(int x, int y, int sublattice, Pauli p);
    void add_errors();
    virtual void ca_correction(int threshold, int sublattice);
};

#endif /* ca_sim_hpp */
