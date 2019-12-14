//
//  Init.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/7/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include "Init.hpp"

using namespace std;
//constructor_default generate sc strucure
Init::Init(int Ncube_, double L_, double T0_, double M_):
    Ncube(Ncube_),
    N(pow(Ncube,3)),
    L(L_),
    T0(T0_),
    M(M_)

{
    Init::alloc_mem(N,dim);
    Init::gen_sc(N,L);
    Init:: gen_ran_vel(N,T0,M,1);
    
}
//constructor specific
Init::Init(string mode, int Ncube_, double L_, double T0_, double M_):
    Ncube(Ncube_),
    L(L_),
    T0(T0_),
    M(M_)
{
    if (mode=="sc" || mode=="liquid"){
        N = pow(Ncube, 3);
        Init::alloc_mem(N,dim);
        Init:: gen_ran_vel(N,T0,M,1);
        
        if (mode =="sc") {
            Init::gen_sc(N,L);
        }
        else Init::gen_rand(N,L);
    }
    
    else if (mode=="fcc"){
        N = pow(Ncube,3) * 4;
        Init::alloc_mem(N,dim);
        Init::gen_fcc(N,L);
        Init:: gen_ran_vel(N,T0,M,1);
    }
    
   
}
//constructor from file
Init::Init(string mode, string file_name, int desired_frame, double T0_, double M_):
    fileName(file_name),
    T0(T0_),
    M(M_)
{
    ifstream ip_file;
    ip_file.open(fileName);
    string ip_word;
    ip_file >> ip_word;
    N = stoi(ip_word);
    L = pow(N,1/3);
    ip_file.close();
    ip_file.clear();
    Init::alloc_mem(N,dim);
    Init::from_file(fileName, desired_frame);
    Init:: gen_ran_vel(N,T0,M,1);
}

//destructor
Init::~Init(){
    
}
//geters
double ** Init::getPosition(){
    return Init::position;
}

double ** Init::getVelocity(){
    return Init::velocity;
}
int Init::getN(){
    return Init::N;
}
double Init::getT(){
    return Init::T0;
}
double Init::getL(){
    return Init::L;
}
int Init::getDim(){
    return Init::dim;
}
int Init::getNcube(){
    return Ncube;
}
//*end geters*//

//memory allocation for position and velocity
void Init::alloc_mem(int N_, int dim_=3){
    position = new double *[N_];
    velocity = new double *[N_];
#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        position[i] = new double[dim_]; // 3-dim
        velocity[i] = new double [dim_]; // 3-dim
    }
#pragma omp critical
    {
        
    }
}

// geneartion routines
void Init::gen_sc(int N_, double L_){
    double rs = L_ / Ncube;
    double roffset = L_ / 2 - rs / 2;
    int n = 0;
#pragma omp paralle for
    for (int x = 0; x < Ncube; ++x) {
        for (int y = 0; y < Ncube; ++y) {
            for (int z = 0; z < Ncube; ++z) {
                if (n < N_){
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
    double rs = L_ / Ncube;
    double roffset = L_ / 2 - rs / 2;
    int n = 0;
#pragma omp paralle for
    for (int x = 0; x < Ncube; ++x) {
        for (int y = 0; y < Ncube; ++y) {
            for (int z = 0; z < Ncube; ++z) {
                if (n < N_){
                    position[n][0] = rs * x - roffset;
                    position[n][1] = rs * y - roffset;
                    position[n][2] = rs * z - roffset;
                }
                n++;
                if (n < N_){
                    position[n][0] = rs * (x + 0.5) - roffset;
                    position[n][1] = rs * (y + 0.5) - roffset;
                    position[n][2] = rs * z - roffset;

                }
                n++;
                if (n < N_){
                    position[n][0] = rs * (x + 0.5) - roffset;
                    position[n][1] = rs * y - roffset;
                    position[n][2] = rs * (z + 0.5) - roffset;

                }
                n++;
                if (n < N_){
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

void Init::gen_rand(int N_, double L_){
    double rs = L_ / Ncube;
    double roffset = L / 2 - rs / 2;
    int n = 0;
    while (n < N_) {
        position[n][0] = (rand()/double(RAND_MAX)) * L_ - roffset;
        position[n][1] = (rand()/double(RAND_MAX)) * L_ - roffset;
        position[n][2] = (rand()/double(RAND_MAX)) * L_ - roffset;
        n++;
        
    }
    
}

void Init::gen_ran_vel(int N_, double T0_, double M_, unsigned int seed = 1){
    //int dim = 3;
    srand(seed);
    double sumV [3] = {0,0,0};
#pragma omp parallel for
    for (int i_atom =0; i_atom < N_; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            velocity[i_atom][i_dim] = (rand()/double(RAND_MAX)) - 0.5;
            sumV [i_dim] += velocity[i_atom][i_dim];
        }
    }
    #pragma omp critical
    {
        
    }
    double KE = 0;
#pragma omp parallel for
    for (int i_atom =0; i_atom < N_; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            velocity[i_atom][i_dim] -= sumV [i_dim];
            KE += velocity[i_atom][i_dim] * velocity[i_atom][i_dim];
        }
    }
    #pragma omp critical
    {
        
    }
    double vscale = pow(((dim * N_ * T0_)/(M_ * KE)), 0.5);
#pragma omp parallel for
    for (int i_atom =0; i_atom < N_; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            velocity[i_atom][i_dim] *= vscale;
        }
    }
    #pragma omp critical
    {

    }
    
    
}


void Init::from_file(string filename, int desired_frame){
    ifstream ip_file;
    ip_file.open(filename);
    string ip_word;
    int frame = 0;
    while (ip_file >> ip_word) {
        getline(ip_file,ip_word);
        getline(ip_file,ip_word);
        for (int i_atom = 0; i_atom < N; ++i_atom) {
            ip_file >> ip_word; // skip Ar
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                ip_file >> ip_word;
                position[i_atom][i_dim] = stod(ip_word);
            }
        }
        if (frame == desired_frame ) {
            break;
        }
        ++frame;
        
    }
    
}
