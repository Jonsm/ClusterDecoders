//
//  sierpinski_sublattice.cpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

#include <iostream>
#include "sierpinski_sublattice.hpp"

using namespace std;

/*
 * Helper class for exact decoder. Each sublattice corresponds to a bipartition of the honeycomb
 * lattice.
 */

sierpinski_sublattice::sierpinski_sublattice(int w, int h) :
w(w),
h(h),
elimination(w),
syndromeMoves(boost::extents[w][h]),
correction(boost::extents[w][h]),
syndromes(boost::extents[w][h]),
error(boost::extents[w][h])
{
    syndromeOffsets  = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    
    getLogicals();
    getEvolutions();
}

//apply Z flip
void sierpinski_sublattice::flip(int x, int y) {
    for (pair<int,int> offset : syndromeOffsets) {
        int shiftX = (x + offset.first + w) % w;
        int shiftY = (y + offset.second + h) % h;
        syndromes[shiftX][shiftY] ^= 1;
    }
    error[x][y] ^= 1;
}

//reset between runs
void sierpinski_sublattice::clear() {
    fill(syndromes.data(), syndromes.data() + syndromes.num_elements(), 0);
    fill(error.data(), error.data() + error.num_elements(), 0);
}

//returns if x,y has been flipped
int sierpinski_sublattice::get(int x, int y) {
    return error[x][y];
}

//initialize logical operators
void sierpinski_sublattice::getLogicals() {
    logicals.push_back(boost::multi_array<int, 2>(boost::extents[w][h]));
    
    for (int i = 0; i < 3; i++) {
        logicals.push_back(boost::multi_array<int, 2>(boost::extents[w][h]));
        boost::multi_array<int, 2>& logical = logicals.back();
        
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                if ((x + i + y) % 3 != 0) {
                    logical[x][y] = 1;
                }
            }
        }
    }
}

//initialize the evolutions of the CA rule. Evolve a single flip on the first row
//to all other rows of the lattice.
void sierpinski_sublattice::getEvolutions() {
    vector<int> r1(w);
    vector<int> r2(w);
    
    vector<int>* currentRow = &r1;
    vector<int>* prevRow = &r2;
    (*currentRow)[0]=1;
    
    for (int y = 0; y < h; y++) {
        vector<int>* tmp = currentRow;
        currentRow = prevRow;
        prevRow = tmp;
        fill(currentRow->begin(), currentRow->end(), 0);
        
        for (int x = 0; x < w; x++) {
            (*currentRow)[x] ^= (*prevRow)[x];
            (*currentRow)[(x-1+w) % w] ^= (*prevRow)[x];
        }
    }
    (*currentRow)[0]^=1;
    
    for (int ofs = 0; ofs < w; ofs++) {
        evolutions.push_back(vector<int>(w));
        for (int x = 0; x < w; x++) {
            evolutions.back()[x] = (*currentRow)[(x - ofs + w) % w];
        }
    }
}

//move all of the error syndromes to one line by applying errors.
void sierpinski_sublattice::moveToLine(boost::multi_array<int, 2> &syndromeMoves) {
    for (int y = h - 1; y > 0; y--) {
        for (int x = 0; x < w; x++) {
            if (syndromeMoves[x][y] == 1) {
                syndromeMoves[x][y-1] ^= 1;
                syndromeMoves[(x-1 + w) % w][y-1] ^= 1;
            }
        }
    }
}

//find the error pattern on the top line given the syndrome information
void sierpinski_sublattice::getTopLine(vector<int>& topLine, boost::multi_array<int, 2> &syndromeMoves) {
    vector<int> topLineSyndrome(w);
    for (int x = 0; x < w; x++) {
        topLineSyndrome[x] = syndromeMoves[x][0];
    }
    if (!elimination.checkAdd(topLineSyndrome.begin(), topLineSyndrome.end())) {
        return;
    }
    
    vector<vector<int>>::iterator evol = evolutions.begin();
    while(elimination.checkAdd(evol->begin(), evol->end())) {
        evol++;
    }
    
    list<int> linearCombination;
    elimination.getLinearCombination(linearCombination);
    for (int i : linearCombination) {
        if (i == 0) {
            continue;
        }
        topLine[i - 1] = 1;
    }
}

//find the most likely correction by checking all possible corrections
void sierpinski_sublattice::mostLikelyCorrection(std::vector<int> &topLine, boost::multi_array<int, 2> &syndromeMoves) {
    for (int x = 0; x < w; x++) {
        correction[x][0] = topLine[x];
    }
    
    for (int y = h - 1; y > 0; y--) {
        for (int x = 0; x < w; x++) {
            correction[x][y] ^= correction[x][(y+1) % h];
            correction[(x-1+w)%w][y] ^= correction[x][(y+1) % h];
        }
    }
    
    for (int x = 0; x < w; x++) {
        for (int y = 1; y < h; y++) {
            correction[x][y] ^= syndromeMoves[x][y];
        }
    }
    
    vector<int> flipWeights(4);
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            for (int i = 0; i < 4; i++) {
                if ((logicals[i][x][y] ^ correction[x][y]) == 1) {
                    flipWeights[i]++;
                }
            }
        }
    }
    
    int minLogical = int(min_element(flipWeights.begin(), flipWeights.end()) - flipWeights.begin());
    if (minLogical != 0) {
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                correction[x][y] ^= logicals[minLogical][x][y];
            }
        }
    }
}

//correct the error
void sierpinski_sublattice::getCorrection() {
    fill(correction.origin(), correction.origin()+correction.num_elements(), 0);
    elimination.clear();
    
    copy(syndromes.data(), syndromes.data() + syndromes.num_elements(), syndromeMoves.data());
    moveToLine(syndromeMoves);
    vector<int> topLine(w);
    getTopLine(topLine, syndromeMoves);
    mostLikelyCorrection(topLine, syndromeMoves);
}

//check if the correction + error is logically trivial
bool sierpinski_sublattice::checkCorrection() {
    getCorrection();
    
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            if (correction[x][y] != error[x][y]) {
                return false;
            }
        }
    }
    return true;
}

//print error
void sierpinski_sublattice::printError() {
    int count = 0;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            count += error[x][y];
        }
    }
    cout << count << "###########" << endl;
    if (count > 0) {
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                cout << error[x][y] << " ";
            }
            cout << "\n";
        }
    }
}

//print violated stabilizers
void sierpinski_sublattice::printSyndrome() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            pair<int, int> coord (x,y);
        }
        cout << "\n";
    }
}

//print correction operator
void sierpinski_sublattice::printCorrection() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << correction[x][y] << " ";
        }
        cout << "\n";
    }
}

//get correction at x,y
int sierpinski_sublattice::correctionAt(int x, int y) {
    return correction[x][y];
}
