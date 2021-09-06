//
//  random_discrete_array.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/22/21.
//

#include "random_discrete_array.hpp"

using namespace std;

//probs = array of weights an item can have. l = number of items. Items start
//assigned to probs[0]
random_discrete_array::random_discrete_array(vector<float> probs, int l) :
n_per_prob(probs.size()),
index_to_prob(l),
index_to_index(l),
dis(0.0, 1.0)
{
    n_probs = (int)probs.size();
    this->probs = probs;
    engine.seed(random_device{}());
    
    for (int i = 0; i < n_probs; i++) {
        prob_to_index.push_back(vector<int>(l));
    }
    for (int i = 0; i < l; i++) {
        prob_to_index[0][i] = i;
    }
    n_per_prob[0] = l;
    for (int i = 0; i < l; i++) {
        index_to_index[i] = i;
    }
}

//change the weight assigned to item i to weight at probs[prob_index]
void random_discrete_array::assign(int i, int prob_index) {
    int former_prob_index = index_to_prob[i];
    int former_index = index_to_index[i];
    int last_at_former = prob_to_index[former_prob_index][n_per_prob[former_prob_index] - 1];
    
    n_per_prob[former_prob_index]--;
    prob_to_index[former_prob_index][former_index] = last_at_former;
    index_to_index[last_at_former] = former_index;
    
    n_per_prob[prob_index]++;
    prob_to_index[prob_index][n_per_prob[prob_index] - 1] = i;
    index_to_prob[i] = prob_index;
    index_to_index[i] =  n_per_prob[prob_index] - 1;
}

//get a random item with probability proportional to its weight
int random_discrete_array::get_random() {
    float R = 0;
    for (int i = 0; i < probs.size(); i++) {
        R += probs[i] * n_per_prob[i];
    }
    float r = dis(engine) * R;

    float R_cumulative = probs[0] * n_per_prob[0];
    int prob_index = 0;
    while (R_cumulative < r) {
        prob_index++;
        R_cumulative += probs[prob_index] * n_per_prob[prob_index];
    }
    
    R_cumulative -= probs[prob_index] * n_per_prob[prob_index];
    float pos = (r - R_cumulative) / (probs[prob_index] * n_per_prob[prob_index]);
    int index = int(pos * n_per_prob[prob_index]);
    if (index == n_per_prob[prob_index]) {
        index--;
    }
    return prob_to_index[prob_index][index];
}

//return the sum of all weights
float random_discrete_array::total_prob() {
    float R = 0;
    for (int i = 0; i < probs.size(); i++) {
        R += probs[i] * n_per_prob[i];
    }
    return R;
}
