//
//  thermal_finite_stabilizers.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 11/1/21.
//

//DEPRECATED

#ifndef thermal_finite_stabilizers_hpp
#define thermal_finite_stabilizers_hpp

#include <vector>
#include <random>
#include "boost/multi_array.hpp"
#include "random_discrete_array.hpp"

class thermal_finite_stabilizers {
public:
    thermal_finite_stabilizers(int w, int h, float T, char stabilizer);
    
    std::pair<int,int> rand_site_thermal();
    void make_flip(int x, int y, bool z);
    float get_R();
    void debug();
private:
    float R;
    const float T;
    const int w;
    const int h;
    const char stabilizer;
    boost::multi_array<int, 2> syndromes;
    boost::multi_array<int, 2> gamma_grid;
    std::vector<std::pair<int,int>> offsets_Z;
    std::vector<std::pair<int,int>> offsets_not_Z;
    std::vector<float> unique_gammas;
    random_discrete_array rd_array;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    
    static std::vector<float> get_unique_gammas(float T);
};

#endif /* thermal_finite_stabilizers_hpp */
