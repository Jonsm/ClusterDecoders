//
//  exact_decoder.hpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

/*
 * Exact decoder, with 50% threshold at infinite bias. Does not work at finite bias. Only works on
 * system sizes where the GSD of the Newman-Moore model is 4. Works by computing all possible corrections,
 * then returning the correction with the lowest weight.
 */

#ifndef exact_decoder_hpp
#define exact_decoder_hpp

#include "paulis.hpp"
#include "sierpinski_sublattice.hpp"

class exact_decoder {
public:
    exact_decoder(int w, int h, std::vector<float> bias, bool give_up=false, bool reduce_weight=false);
    void flip(int x, int y, int sublattice, Pauli p);
    bool make_correction();
    int num_errors();
    
    //debug functions
    void print_error();
    void print_correction();
    void print_syndrome();
    
    //reset errors
    void clear();
private:
    const int w;
    const int h;
    sierpinski_sublattice lattice_A;
    sierpinski_sublattice lattice_B;
};

#endif /* exact_decoder_hpp */
