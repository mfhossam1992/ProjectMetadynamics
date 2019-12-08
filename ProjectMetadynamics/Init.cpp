//
//  Init.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/7/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include "Init.hpp"

using namespace std;
//constructor_default
Init::Init(int Ncube_, double L_):
    Ncube(Ncube_),
    N(pow(Ncube,3)),
    L(L_)

{
    Init::alloc_mem(N);
    Init::gen_sc(N,L);
    
}
//constructor specific
Init::Init(string mode, int Ncube_, double L_):
    Ncube(Ncube_),
    L(L_)
{
    if (mode=="sc" || mode=="liquid"){
        N = pow(Ncube, 3);
        Init::alloc_mem(N);
        
        if (mode =="sc") {
            Init::gen_sc(N,L);
        }
        //else Init::gen_rand(N,L);
    }
    
    else if (mode=="fcc"){
        N = pow(Ncube,3) * 4;
        Init::alloc_mem(N);
        Init::gen_fcc(N,L);
    }
   
}
//destructor
Init::~Init(){
    
}
//geters
double ** Init::getPosition(){
    return Init::position;
}

int Init::getN(){
    return Init::N;
}
//*end geters*//

//memory allocation for position and velocity
void Init::alloc_mem(int N_){
    position = new double *[N_];
    velocity = new double *[N_];
#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        position[i] = new double[3];
        velocity[i] = new double [3];
    }
#pragma omp critical
    {
        
    }
}

// geneartion routines
void Init::gen_sc(int N_, double L_){
    double rs = L / Ncube;
    double roffset = L / 2 - rs / 2;
    int n = 0;
#pragma omp paralle for
    for (int x = 0; x < Ncube; ++x) {
        for (int y = 0; y < Ncube; ++y) {
            for (int z = 0; z < Ncube; ++z) {
                if (n < N){
                    position[n][0] = rs * x - roffset;
                    position[n][1] = rs * y - roffset;
                    position[n][2] = rs * z - roffset;
                }
                n++;
            }
        }
    }
#pragma omp critical
    {
        
    }
    
}

void Init::gen_fcc(int N_, double L_){
    double rs = L / Ncube;
    double roffset = L / 2 - rs / 2;
    int n = 0;
#pragma omp paralle for
    for (int x = 0; x < Ncube; ++x) {
        for (int y = 0; y < Ncube; ++y) {
            for (int z = 0; z < Ncube; ++z) {
                if (n < N){
                    position[n][0] = rs * x - roffset;
                    position[n][1] = rs * y - roffset;
                    position[n][2] = rs * z - roffset;
                }
                n++;
                if (n < N){
                    position[n][0] = rs * (x + 0.5) - roffset;
                    position[n][1] = rs * (y + 0.5) - roffset;
                    position[n][2] = rs * z - roffset;

                }
                n++;
                if (n < N){
                    position[n][0] = rs * (x + 0.5) - roffset;
                    position[n][1] = rs * y - roffset;
                    position[n][2] = rs * (z + 0.5) - roffset;

                }
                n++;
                if (n < N){
                    position[n][0] = rs * x  - roffset;
                    position[n][1] = rs * (y + 0.5) - roffset;
                    position[n][2] = rs * (z + 0.5) - roffset;

                }
                n++;


            }
        }
    }
#pragma omp critical
    {
        
    }
    
}
