//
//  exact_decoder.cpp
//  ClusterDecoders
//
//  Created by Jon on 1/4/22.
//

#include <iostream>
#include "exact_decoder.hpp"

using namespace std;

exact_decoder::exact_decoder(int w, int h, vector<float> bias, bool give_up, bool reduce_weight) :
w(w),
h(h),
lattice_A(w, h),
lattice_B(w, h)
{
    
}

//flip one qubit, should be Z or it fails
void exact_decoder::flip(int x, int y, int sublattice, Pauli p) {
    if (p != Z) {
        cout << "PARTIAL BIAS USED WITH EXACT DECODER" << endl;
        return;
    }
    
    if (sublattice == 1) {
        lattice_A.flip(x, y);
    } else {
        lattice_B.flip(w-x-1, h-y-1);
    }
}

int exact_decoder::num_errors() {
    return -1; //implement
}

//make correction
bool exact_decoder::make_correction() {
    return lattice_A.checkCorrection() && lattice_B.checkCorrection();
}

//reset all errors
void exact_decoder::clear() {
    lattice_A.clear();
    lattice_B.clear();
}

//print errors
void exact_decoder::print_error() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (lattice_A.error[x][y]) {
                cout << "Z";
            } else {
                cout << "-";
            }
            if (lattice_B.error[w-x-1][h-y-1]) {
                cout << "Z";
            } else {
                cout << "-";
            }
            cout << " ";
        }
        cout << endl;
    }
}

//print violated stabilizers
void exact_decoder::print_syndrome() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            pair<int,int> coord_A(x,y);
            cout << lattice_A.syndromes[coord_A.first][coord_A.second];
            pair<int,int>coord_B(w-x-1, (2*h - y - 2) % h);
            cout << lattice_A.syndromes[coord_B.first][coord_B.second];
            cout << " ";
        }
        cout << endl;
    }
}

//print correction operator
void exact_decoder::print_correction() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            if (lattice_A.correction[x][y]) {
                cout << "Z";
            } else {
                cout << "-";
            }
            if (lattice_B.correction[w-x-1][h-y-1]) {
                cout << "Z";
            } else {
                cout << "-";
            }
            cout << " ";
        }
        cout << endl;
    }
}
