//
//  ClusterDecoder.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 8/4/21.
//

#include "cluster_decoder.hpp"
#include <iostream>

using namespace std;

cluster_decoder::cluster_decoder(int w, int h, vector<float> bias, bool give_up, bool reduce_weight) :
w(w),
h(h),
lattice_A(w,h,bias, 'A', give_up, reduce_weight),
lattice_B(w,h,bias, 'B', give_up, reduce_weight)
{
    this->bias = bias;
}

void cluster_decoder::flip(int x, int y, int sublattice, Pauli p) {
    if (sublattice == 1) {
        switch (p) {
            case X:
                lattice_B.flip(x, y, false);
                break;
            case Y:
                lattice_B.flip(x, y, false);
                lattice_A.flip(x, y, true);
                break;
            case Z:
                lattice_A.flip(x, y, true);
                break;
        }
    } else {
        switch (p) {
            case X:
                lattice_A.flip(x, y, false);
                lattice_B.flip(x, y, true);
                break;
            case Y:
                lattice_A.flip(x, y, false);
                break;
            case Z:
                lattice_B.flip(x, y, true);
                break;
        }
    }
}

bool cluster_decoder::make_correction() {
    return lattice_A.make_correction() && lattice_B.make_correction();
}

bool cluster_decoder::check_correction() {
    return lattice_A.check_correction() && lattice_B.check_correction();
}

void cluster_decoder::clear() {
    lattice_A.clear();
    lattice_B.clear();
}

void cluster_decoder::print_helper(boost::multi_array<int, 2>& z_A, boost::multi_array<int, 2>& not_Z_A,boost::multi_array<int, 2>& z_B, boost::multi_array<int, 2>& not_Z_B) {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int z1 = z_A[x][y];
            int x1 = not_Z_B[x][y];
            int z2 = z_B[x][y];
            int y2 = not_Z_A[x][y];
            
            if (z1 && x1) {
                cout << "Y";
            } else if (z1 && !x1) {
                cout << "Z";
            } else if (!z1 && x1) {
                cout << "X";
            } else {
                cout << "-";
            }
            
            if (z2 && y2) {
                cout << "X";
            } else if (z2 && !y2) {
                cout << "Z";
            } else if (!z2 && y2) {
                cout << "Y";
            } else {
                cout << "-";
            }
            
            cout << " ";
        }
        cout << endl;
    }
}

void cluster_decoder::print_error() {
    print_helper(lattice_A.errors_Z, lattice_A.errors_not_Z, lattice_B.errors_Z, lattice_B.errors_not_Z);
}

void cluster_decoder::print_syndrome() {
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << lattice_A.syndromes[x][y] << lattice_B.syndromes[x][y] << " ";
        }
        cout << endl;
    }
}

void cluster_decoder::print_correction() {
    print_helper(lattice_A.correction_Z, lattice_A.correction_not_Z, lattice_B.correction_Z, lattice_B.correction_not_Z);
}

void cluster_decoder::print_error_after_correction() {
    boost::multi_array<int, 2> total_A_Z(boost::extents[w][h]);
    boost::multi_array<int, 2> total_A_not_Z(boost::extents[w][h]);
    boost::multi_array<int, 2> total_B_Z(boost::extents[w][h]);
    boost::multi_array<int, 2> total_B_not_Z(boost::extents[w][h]);
    for (int x = 0; x < w; x++) {
        for (int y = 0; y < h; y++) {
            total_A_Z[x][y] = lattice_A.correction_Z[x][y] ^ lattice_A.errors_Z[x][y];
            total_A_not_Z[x][y] = lattice_A.correction_not_Z[x][y] ^ lattice_A.errors_not_Z[x][y];
            total_B_Z[x][y] = lattice_B.correction_Z[x][y] ^ lattice_B.errors_Z[x][y];
            total_B_not_Z[x][y] = lattice_B.correction_not_Z[x][y] ^ lattice_B.errors_not_Z[x][y];
        }
    }
    
    print_helper(total_A_Z, total_A_not_Z, total_B_Z, total_B_not_Z);
}
