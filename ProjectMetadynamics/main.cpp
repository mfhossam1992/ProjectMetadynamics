//
//  main.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/7/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include <iostream>
#include"Init.hpp"

using namespace std;

// Input Parameters

// ************* System ************* //

// mass
double M = 48.0;
// Number of Particles Per dimension
int Ncube = 4;
// Box Side Length
double L = 1.56 * Ncube;
// Initial Temperature
double T0 = 1;
// System Temperature
double Ta = 2;

// ************* MD ************* //

// cutoff radius
double rc = L/2;
//Anderson Thermostate
bool anderson = true;
double eta = 0.3125;
// Output Atomic Trajcetory filename
string fileName = "out_traj";
// total numern of time steps
int steps = 5000;
// time step size
double h = 0.032;

// ************* MeTaD ************* //

//Gaussian
double meta_w = 0.1; // height
double meta_sigma = 0.1; // width
// Maximum Number of Gaussian
int meta_max = 100;
// frequency
int meta_tau = 50;
// cutoff radius for n.n. in Q6
double meta_rc = 1.2 * pow(2, 1/6); // 1.2 * r_min in LJ


// ** end Input Parameters ** //

int main(int argc, const char * argv[]) {
    
    Init init("fcc",Ncube,L);
    int N = init.getN();
    cout<< N<< endl;
    return 0;
}
