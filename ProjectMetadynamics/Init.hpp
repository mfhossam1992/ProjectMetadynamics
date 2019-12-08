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


#ifdef OMP
#include "omp.h"
#endif

using namespace std;

class Init {
    int Ncube;
    int N ;
    double L = 4;
    double ** position;
    double ** velocity;
    void alloc_mem(int);
    void gen_sc(int,double);
    void gen_fcc(int,double);
    //void gen_rand(int,double);
    //void gen_ran_vel(float);
    
    
    
    
public:
    //constructor(s)
    Init(int, double);
    Init(string, int, double);
    //destructor
    virtual ~Init();
    double ** getPosition();
    double ** getVelocity();
    int getN();
};
#endif /* Init_hpp */
