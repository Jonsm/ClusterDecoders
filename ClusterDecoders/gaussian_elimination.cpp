//
//  gaussian_elimination.cpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

#include "gaussian_elimination.hpp"

using namespace std;

//initialize gaussian elimination on Z2 vectors of length l
gaussian_elimination::gaussian_elimination(int l) : l(l), vecTmp(l) {
    
}

//Add a vector to the current list of vectors. If the vector is linearly dependent, return false.
//begin and end are iterators in the source array, from which to copy. end should be begin+l.
bool gaussian_elimination::checkAdd(vector<int>::iterator begin, vector<int>::iterator end) {
    vecTmp.assign(begin, end);
    lastIndAdded++;
    
    for (pair<int, int> rowOp : rowOps) {
        vecTmp[rowOp.second] = vecTmp[rowOp.first] ^ vecTmp[rowOp.second];
    }
    
    int foundPivot = -1;
    for (int i = 0; i < vecTmp.size(); i++) {
        if (vecTmp[i] == 1 && pivots.find(i) == pivots.end()) {
            pivots.insert(i);
            pivotsList.push_back(i);
            foundPivot = i;
            break;
        }
    }
    
    if (foundPivot != -1) {
        for (int i = 0; i < vecTmp.size(); i++) {
            if (vecTmp[i] == 1 && i != foundPivot) {
                rowOps.push_back(pair<int,int>(foundPivot, i));
            }
        }
        return true;
    } else {
        return false;
    }
}

//return the linear combination of vectors equal to the last vector added. Only works
//if the last call of checkAdd is false.
void gaussian_elimination::getLinearCombination(list<int>& linearCombination) {
    int i = 0;
    for (list<int>::iterator it = pivotsList.begin(); it != pivotsList.end(); it++) {
        int pivot = *it;
        if (vecTmp[pivot] == 1) {
            linearCombination.push_back(i);
        }
        i++;
    }
    
    linearCombination.push_back(lastIndAdded - 1);
}

//reset the gaussian elimination, clear all vectors from the list
void gaussian_elimination::clear() {
    rowOps.clear();
    pivots.clear();
    pivotsList.clear();
    lastIndAdded = 0;
}
