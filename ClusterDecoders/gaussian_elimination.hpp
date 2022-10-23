//
//  gaussian_elimination.hpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

/*
 * Class for Gaussian elimination mod 2. Used to solve for error given syndrome.
 */

#ifndef gaussian_elimination_hpp
#define gaussian_elimination_hpp

#include <vector>
#include <set>
#include <list>

class gaussian_elimination {
private:
    int l;
    int lastIndAdded = 0;
    std::vector<int> vecTmp;
    std::list<std::pair<int,int>> rowOps;
    std::set<int> pivots;
    std::list<int> pivotsList;
    
public:
    gaussian_elimination(int l);
    bool checkAdd(std::vector<int>::iterator begin, std::vector<int>::iterator end);
    void getLinearCombination(std::list<int>& linearCombination);
    void clear();
};

#endif /* gaussian_elimination_hpp */
