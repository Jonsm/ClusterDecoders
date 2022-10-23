//
//  random_ca_sim.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 10/26/21.
//

#ifndef random_ca_sim_hpp
#define random_ca_sim_hpp

#include <vector>
#include <random>
#include "ca_sim.hpp"

class random_ca_sim : public ca_sim {
public:
    random_ca_sim(int w, int h, std::vector<float>& bias);
    
    static void distribution_sim(int avg, int t, std::vector<std::pair<int,int>>& dims, std::vector<float>& ps, std::vector<float>& normalized_bias, std::string filename);
private:
    const float p;
    float beta;
    float dt;
    std::vector<float> fs;
    
    void init_constants();
    void ca_correction(int threshold, int sublattice) override;
};

#endif /* random_ca_sim_hpp */
