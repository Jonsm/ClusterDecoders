//
//  random_discrete_array.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/22/21.
//

//helper data structure for thermal_bias_sim. array that supports
// 1) assigning weights to elements from a predefined list probs
// 2) getting a random element with probability proportional to its weight

#ifndef random_discrete_array_hpp
#define random_discrete_array_hpp

#include <vector>
#include <random>

class random_discrete_array {
public:
    random_discrete_array(std::vector<float> probs, int l);
    void assign(int i, int prob_index);
    int get_random();
    float total_prob();
    int get_prob_index(int i);
    
private:
    int l;
    int n_probs;
    std::vector<float> probs;
    std::vector<std::vector<int>> prob_to_index;
    std::vector<int> n_per_prob;
    std::vector<int> index_to_prob;
    std::vector<int> index_to_index;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
};

#endif /* random_discrete_array_hpp */
