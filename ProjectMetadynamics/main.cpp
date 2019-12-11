//
//  main.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/7/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include <iostream>
#include"Init.hpp"
#include "MD.hpp"

using namespace std;

// Input Parameters

// ************* System ************* //

// mass
double M = 48.0;
// Number of Particles Per dimension
int Ncube = 4 ;
// Box Side Length
//double L = 4;
double L = 1.56 * Ncube;
//double L = 1.49 * Ncube;
// Initial Temperature
double T0 = 0.7;
// System Temperature
double Ta = 0.8;

// ************* MD ************* //

// cutoff radius
double rc = L/2;//
//Anderson Thermostate
bool anderson = true;
double eta = 0.3125;
// Output Atomic Trajcetory filename
string fileName = "thermo.txt";
string trajFileName = "md.xyz";
// total numern of time steps
int steps = 50000;
// time step size
double h = 0.032;//

// ************* MeTaD ************* //

//Gaussian
double meta_w = 0.1; // height
double meta_sigma = 0.1; // width of Q_6
double meta_sigma_2 = 0.1; // width of potential energy
int cv = 2; // change this to be 1 in case of only 1 CV (Q6)
// Maximum Number of Gaussian
int meta_max = 1000;
// frequency
int meta_tau = 50;
// cutoff radius for n.n. in Q6
double meta_rc = 1.2 * pow(2, 1/6); // 1.2 * r_min in LJ
// perform or not mtd
bool mtd = true;

// ** end Input Parameters ** //

int main(int argc, const char * argv[]) {
    
    Init * init = new Init("fcc",Ncube,L,T0,M);
   // int desired_frame_number = 11;
   // Init * init = new Init("file","/Users/hossamfarag/Desktop/Hossam/MSE485 - Fa2019/Project/metadynamics/ProjectMetadynamics/ProjectMetadynamics/Debug/init.xyz",desired_frame_number, T0,M);
    
    auto N = init->getN();
    cout<< "N = " << N << endl;
    auto T = init->getT();
    cout << "T0 = " << T << endl;
    cout << "Ta = " << Ta << endl;
    cout << "BOX LENGTH:"<< L <<endl;
    auto R = init->getPosition();
    cout << "\n Positions\n" << endl;
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        cout << R[i_atom][0] << " , " << R[i_atom][1] << " , " << R[i_atom][2] << endl;

    }
    auto V = init->getVelocity();
    cout << "\n Velocities\n" << endl;
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        cout << V[i_atom][0] << " , " << V[i_atom][1] << " , " << V[i_atom][2] << endl;

    }

    if (mtd == false) {
        MD md(init, anderson, Ta, eta, mtd, rc, meta_rc, h, fileName, steps,trajFileName );
        md.simulate();

    } else {
        if (cv == 2) {
            MD md(init, anderson, Ta, eta, mtd, meta_w, meta_sigma, meta_sigma_2,meta_max, meta_tau,rc, meta_rc, h, fileName, steps,trajFileName );
            md.simulate();
        } else {
            MD md(init, anderson, Ta, eta, mtd, meta_w, meta_sigma, meta_max, meta_tau,rc, meta_rc, h, fileName, steps,trajFileName );
            md.simulate();
        }

    }

    return 0;
}
