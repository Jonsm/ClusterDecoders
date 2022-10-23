//
//  thermal_finite_stabilizers.cpp
//  ClusterDecoders
//
//  Created by Jon San Miguel on 11/1/21.
//

#include <iostream>
#include "thermal_finite_stabilizers.hpp"

using namespace std;

vector<float> thermal_finite_stabilizers::get_unique_gammas(float T) {
    std::vector<float> unique_gammas;
    double beta = 1.0 / T;
    for (double omega = -3; omega <= 3; omega += 2) {
        double gamma = omega / (1.0 - exp(-beta * omega));
        unique_gammas.push_back(gamma);
//        cout << "G:" << gamma << endl;
    }
    return unique_gammas;
}

thermal_finite_stabilizers::thermal_finite_stabilizers(int w, int h, float T, char stabilizer) :
w(w),
h(h),
T(T),
stabilizer(stabilizer),
syndromes(boost::extents[w][h]),
gamma_grid(boost::extents[w][h]),
rd_array(get_unique_gammas(T), w*h),
dis(0.0,1.0)
{
    engine.seed(std::random_device{}());
    vector<pair<int,int>> s1 = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(-1,-1)};
    vector<pair<int,int>> s2 = vector<pair<int,int>> {pair<int,int>(0,0),pair<int,int>(0,-1),pair<int,int>(1,0)};
    if (stabilizer == 'A') {
        offsets_Z = s1;
        offsets_not_Z = s2;
    } else {
        offsets_Z = s2;
        offsets_not_Z = s1;
    }
    
    R = rd_array.total_prob();
}

void thermal_finite_stabilizers::debug() {
    cout << stabilizer << ":::::::::" << endl;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << syndromes[x][y];
        }
        cout << endl;
    }
    
    int counts[] {0,0,0,0};
    
    cout << "Gammagrid:" << endl;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            cout << gamma_grid[x][y];
            counts[gamma_grid[x][y]] += 1;
        }
        cout << endl;
    }
    
    for (int i = 0; i < 4; i++) {
        cout << counts[i] << " ";
    }
    cout << endl;
    
    cout << rd_array.total_prob() << "=" << R << endl;
}

float thermal_finite_stabilizers::get_R() {
    return R;
}

pair<int,int> thermal_finite_stabilizers::rand_site_thermal() {
    int ind = rd_array.get_random();
    int x = ind % w;
    int y = ind / w;
    return pair<int, int>(x,y);
}

void thermal_finite_stabilizers::make_flip(int x, int y, bool z) {
    vector<pair<int,int>>* offsets = &offsets_Z;
    if (!z) {
        offsets = &offsets_not_Z;
    }
    
    for (const pair<int, int>& coord2 : *offsets) {
        int syndrome_before = syndromes[(x + coord2.first + w) % w][(y + coord2.second + h) % h];
        int x2 = (x + coord2.first + w) % w;
        int y2 = (y + coord2.second + h) % h;
        syndromes[x2][y2] ^= 1;
        
        for (const pair<int, int>& coord3 : offsets_Z) {
            int x3 = (x2 - coord3.first + w) % w;
            int y3 = (y2 - coord3.second + h) % h;
            int gg_before = gamma_grid[x3][y3];
            pair<int, int> coord4(x3,y3);
            
            int gg_new = gg_before + (1 - 2 * syndrome_before);
            gamma_grid[x3][y3] = gg_new;
            int ind = y3 * w + x3;
            rd_array.assign(ind, gg_new);
        }
    }
    R = rd_array.total_prob();
}
