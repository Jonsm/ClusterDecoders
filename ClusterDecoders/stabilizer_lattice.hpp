//
//  StabilizerLattice.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

//helper class for cluster_decoder. Each lattice decodes separately. There
//are two stabilizers A and B, representing (going around hexagon) XZXZXZ and
//ZYZYZY. This is equivalent under local unitary to decoding Z and X errors in CSS code.

#ifndef stabilizer_lattice_hpp
#define stabilizer_lattice_hpp

#include <vector>
#include "clustering_helper.hpp"
#include "bias_correction_helper.hpp"

class stabilizer_lattice {
public:
    boost::multi_array<int, 2> errors_Z;
    boost::multi_array<int, 2> errors_not_Z;
    boost::multi_array<int, 2> syndromes;
    boost::multi_array<int, 2> correction_Z;
    boost::multi_array<int, 2> correction_not_Z;
    
    stabilizer_lattice(int w, int h, std::vector<float> bias, char stabilizer, bool give_up, bool reduce_weight);
    void flip(int x, int y, bool z);
    bool make_correction();
    bool check_correction();
    void clear();
private:
    bool reduce_weight;
    bool give_up;
    int w;
    int h;
    std::vector<float> bias;
    char stabilizer;
    std::vector<std::pair<int, int>> offsets_Z;
    std::vector<std::pair<int, int>> offsets_not_Z;
    clustering_helper clustering;
    boost::multi_array<int, 2> syndromes_tmp;
    bias_correction_helper bias_correction;
    
    void get_crossover_n();
    void attempt_correction(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    bool no_bias_correctable(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    bool correctable(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    bool bias_correctable(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    void string_clean(int x_src, int y_src, int x_dest, int y_dest);
    std::pair<int,int> string_clean_x(int x_src, int y_src, int x_dest);
    void string_clean_y(int x_src, int y_src, int x_dest, int y_dest);
    void no_bias_correction(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    bool check_stabilizers();
    bool check_logicals();
    void reduce_stabilizer_weight();
};

#endif /* StabilizerLattice_hpp */
