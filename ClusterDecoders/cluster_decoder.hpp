//
//  ClusterDecoder.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

//decoder based on RG decoder of Bravyi and Haah. Attempts to correct Pauli-Z errors using
//fractal operators, then corrects remaining errors with stringlike operators.

#ifndef cluster_decoder_hpp
#define cluster_decoder_hpp

#include <vector>
#include "boost/multi_array.hpp"
#include "paulis.hpp"
#include "stabilizer_lattice.hpp"

class cluster_decoder {
public:
    //constructor. give_up = whether decoder gives up after clusters are bigger than L/2. reduce_weight =
    //try to reduce weight of final corrections by multiplying with stabilizers. bias parameter currently
    //unused.
    cluster_decoder(int w, int h, std::vector<float> bias, bool give_up=true, bool reduce_weight=false);
    
    //apply Pauli P to qubit.
    void flip(int x, int y, int sublattice, Pauli p);
    
    //attempt correction. if give_up is true, returns false if the clusters are larger than L/2.
    bool make_correction();
    
    //return true if the correction + error is trivial.
    bool check_correction();
    
    //debug functions
    void print_error();
    void print_correction();
    void print_syndrome();
    void print_error_after_correction();
    
    //reset errors
    void clear();
private:
    int w;
    int h;
    std::vector<float> bias;
    stabilizer_lattice lattice_A;
    stabilizer_lattice lattice_B;
    
    void print_helper(boost::multi_array<int, 2>& z_A, boost::multi_array<int, 2>& not_Z_A,boost::multi_array<int, 2>& z_B, boost::multi_array<int, 2>& not_Z_B);
};

#endif /* ClusterDecoder_hpp */
