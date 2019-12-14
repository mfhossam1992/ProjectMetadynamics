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
int Ncube = 3 ;
// Box Side Length
//double L = 4;
double L = 1.54 * Ncube;
//double L = 1.4938 * Ncube; // corresponding to rho = 1.2
// Initial Temperature
double T0 = 0.7;
// System Temperature
double Ta = 0.7;

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
int steps = 1000000;
// time step size
double h = 0.032;//

// ************* MeTaD ************* //

//Gaussian
double meta_w = 1; // height
double meta_sigma = 0.025; // width of Q_6
double meta_sigma_2 = 15; // width of potential energy // put it to zero if not want to activate it
// Maximum Number of Gaussian
int meta_max = 10000;
// frequency
int meta_tau = 100;
// cutoff radius for n.n. in Q6
double meta_rc = 1.2 * pow(2, 1/6); // 1.2 * r_min in LJ
// perform or not mtd
bool mtd = true;

// ** end Input Parameters ** //

int main(int argc, const char * argv[]) {
    
    // Initiate a box of Ar atoms with the desired number of atoms, structure and Temperature
    Init * init = new Init("fcc",Ncube,L,T0,M);
    //Printing to screen the Parameters of the object init I created
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
    //END of Printing Section
    //Equilibration Section
    double equil_Temp = 0.8;
    double equil_steps = 20000;
    string equil_thermo_file = "equil_thermo.txt";
    string equil_traj_file = "equil_md.xyz";
    MD * equilib = new MD(init, anderson, equil_Temp, eta, false, rc, meta_rc, h, equil_thermo_file, equil_steps, equil_traj_file );
    //Printing to screen the Parameters of the object equilib I created
    N = equilib->getN();
    cout<< "N = " << N << endl;
    double L_ = equilib->getL();
    cout << "BOX LENGTH:"<< L_ <<endl;
    R = equilib->getPosition();
    cout << "\n Positions\n" << endl;
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        cout << R[i_atom][0] << " , " << R[i_atom][1] << " , " << R[i_atom][2] << endl;

    }
    V = equilib->getVelocity();
    cout << "\n Velocities\n" << endl;
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        cout << V[i_atom][0] << " , " << V[i_atom][1] << " , " << V[i_atom][2] << endl;

    }
    //END of Printing Section
    equilib->simulate();
    //End of Equilibration section
    if (mtd == false) {
        MD md(equilib, anderson, Ta, eta, mtd, rc, meta_rc, h, fileName, steps,trajFileName );
        md.simulate();

    }
    else {
        MD md(equilib, anderson, Ta, eta, mtd, meta_w, meta_sigma, meta_sigma_2,meta_max, meta_tau,rc, meta_rc, h, fileName, steps,trajFileName );
        md.simulate();

    }

    return 0;
}
