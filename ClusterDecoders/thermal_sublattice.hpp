//
//  thermal_sublattice.hpp
//  ClusterDecoders
//
//  Created by Jon on 12/1/21.
//

/*
 * Helper for thermal_bias_sim_full. Runs BKL algorithm for a single sublattice.
 */

#ifndef thermal_sublattice_hpp
#define thermal_sublattice_hpp

#include <vector>
#include "random_discrete_array.hpp"
#include "paulis.hpp"

class thermal_sublattice {
public:
    thermal_sublattice(int w, int h, int sublattice, std::vector<float>& bias, bool no_correct=false, bool no_correct_bias=false);
    float R_total();
    Pauli flip_type();
    std::pair<int,int> flip_coord(Pauli pauli);
    void update_probs(int x, int y, char stabilizer, int old_st);
    void debug();
    
private:
    const int w;
    const int h;
    const int sublattice;
    random_discrete_array x_rdarray;
    random_discrete_array y_rdarray;
    random_discrete_array z_rdarray;
    std::mt19937 engine;
    std::uniform_real_distribution<float> dis;
    
    static std::vector<float> get_unique_gammas(std::vector<float>& bias, Pauli pauli, bool no_correct=false, bool no_correct_bias=false);
};

#endif /* thermal_sublattice_hpp */
