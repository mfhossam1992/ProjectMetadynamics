//
//  Init.hpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/7/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#ifndef Init_hpp
#define Init_hpp

#include <stdio.h>
#include <cmath>
#include <string>
#include <random>
#include <fstream>


#ifdef OMP
#include "omp.h"
#endif

using namespace std;

class Init {
    //private members
    string fileName_r; // in case of loading from file (position)
    string fileName_v; // in case of loading from file (velocities)
    int Ncube;
    int N ;
    int dim = 3; // default
    double L = 4; // default = 4
    double T0 = 1; // default 
    double M = 48;
    // trajectories
    double ** position;
    double ** velocity;
    //private methods
    void alloc_mem(int,int);
    void gen_sc(int,double);
    void gen_fcc(int,double);
    void gen_rand(int,double);
    void gen_ran_vel(int, double, double, unsigned int); // N , T0, M, seed
    void from_file_r(string, int); // filename, desired frame number // this function reads positions
    void from_file_v(string, int); //filename, desired frame number // this function reads velocities
    
    
    
    
public:
    //constructor(s)
    Init(int, double, double, double); // Ncube, L, T0, M
    Init(string, int, double, double, double); // mode, Ncube, L, T0, M
    Init(string, string, int, double, double); // mode, filename, desired_frame_number,  T0, M
    Init(string, string, string, int); // mode, filename_r, filename_v, desired_frame_number

    //destructor
    virtual ~Init();
    //getters
    double ** getPosition();
    double ** getVelocity();
    int getN();
    double getT();
    double getL();
    int getDim();
    int getNcube();

};
#endif /* Init_hpp */
