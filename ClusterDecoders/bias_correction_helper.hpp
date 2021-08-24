//
//  bias_correction_helper.hpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/16/21.
//

#ifndef bias_correction_helper_hpp
#define bias_correction_helper_hpp

#include "boost/multi_array.hpp"
#include <vector>

class bias_correction_helper {
public:
    bias_correction_helper(int w, int h, char stabilizer);
    bool check_cluster(std::vector<std::pair<int,int>>& cluster, std::vector<int>& rect);
    void correct_cluster(std::vector<std::pair<int,int>>& cluster, boost::multi_array<int, 2>& correction_Z);
    
private:
    int w;
    int h;
    int rect_w;
    int rect_h;
    int tl_x;
    int tl_y;
    char stabilizer;
    boost::multi_array<int,3> local_gauss_laws;
    boost::multi_array<int,2> corrected_syndromes;
    std::vector<int> current_gl;
    std::vector<int> pivots;
    std::vector<std::pair<int,int>> row_ops;
    
    void get_gauss_laws();
    boost::multi_array<int, 3>::array_view<3>::type get_gauss_law_view(std::vector<std::pair<int,int>>& cluster);
    bool gaussian_elimination(std::vector<std::pair<int,int>>& cluster);
    void apply_correction(boost::multi_array<int, 2>& correction_Z);
    void update_cluster(std::vector<std::pair<int,int>>& cluster);
};

#endif /* bias_correction_helper_hpp */
