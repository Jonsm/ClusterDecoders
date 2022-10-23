//
//  sierpinski_sublattice.hpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

#ifndef sierpinski_sublattice_hpp
#define sierpinski_sublattice_hpp

#include "boost/multi_array.hpp"
#include "gaussian_elimination.hpp"
#include <vector>
#include <set>

class sierpinski_sublattice {
private:
    int w;
    int h;
    std::vector<std::pair<int, int>> syndromeOffsets;
    std::vector<boost::multi_array<int, 2>> logicals;
    std::vector<std::vector<int>> evolutions;
    boost::multi_array<int, 2> syndromeMoves;
    gaussian_elimination elimination;
    
    void getLogicals();
    void getEvolutions();
    void getCorrection();
    void moveToLine(boost::multi_array<int,2>& syndromeMoves);
    void getTopLine(std::vector<int>& topLine, boost::multi_array<int, 2> &syndromeMoves);
    void mostLikelyCorrection(std::vector<int>& topLine, boost::multi_array<int, 2> &syndromeMoves);
    
public:
    boost::multi_array<int, 2> error;
    boost::multi_array<int, 2> correction;
//    std::set<std::pair<int,int>> syndromes;
//    std::vector<std::pair<int,int>> syndromes;
//    int n_syndromes;
    boost::multi_array<int, 2> syndromes;
    
    sierpinski_sublattice(int w, int h);
    void flip(int x, int y);
    int get(int x, int y);
    void clear();
    bool checkCorrection();
    void printError();
    void printCorrection();
    void printSyndrome();
    int correctionAt(int x, int y);
};

#endif /* sierpinski_sublattice_hpp */
