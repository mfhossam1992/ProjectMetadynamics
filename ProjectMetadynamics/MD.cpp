//
//  MD.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/8/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include "MD.hpp"

using namespace std;

// Constructor
MD::MD(Init*& init_, bool anderson_, bool mtd_, double rc_, double metarc_, double h_):
    N(init_->getN()),
    dim(init_->getDim()),
    L(init_->getL()),
    T0(init_->getT()),
    anderson(anderson_),
    mtd(mtd_),
    rc(rc_),
    metarc(metarc_),
    h(h_)

{
    alloc_mem(N, dim);
    alloc_mem_rv(N, dim);
    R = init_->getPosition();
    V = init_->getVelocity();
    alloc_mem3(my_displacement_table_, N, N, dim);
    alloc_mem2(my_distance_table_, N, N);
    my_force_on_ = new double [dim];
    
}
//Destructor
MD::~MD(){
    
}

// memory Allocation Functions

//memory allocation for position and velocity
void MD::alloc_mem(int N_, int dim_=3){
    R = new double *[N_];
    V = new double *[N_];
#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        R[i] = new double[dim_]; // 3-dim
        V[i] = new double [dim_]; // 3-dim
    }
#pragma omp critical
    {
        
    }
}
//memory allocation for new position and new velocity
void MD::alloc_mem_rv(int N_, int dim_=3){
    nR = new double *[N_];
    nV = new double *[N_];
#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        nR[i] = new double[dim_]; // 3-dim
        nV[i] = new double [dim_]; // 3-dim
    }
#pragma omp critical
    {
        
    }
}

//memory allocation for passed double pointer
void MD::alloc_mem2(double ** & dptr, int dim1, int dim2){
    dptr = new double *[dim1];
#pragma omp parallel for
    for (int i =0; i<dim1; ++i) {
        dptr[i] = new double[dim2];
    }
    #pragma omp critical
    {
        
    }
    
}
//memory allocation for passed triple pointer
void MD::alloc_mem3(double *** & dptr, int dim1, int dim2, int dim3){
    dptr = new double **[dim1];
#pragma omp parallel for
    for (int i =0; i<dim1; ++i) {
        dptr[i] = new double * [dim2];
        for (int j = 0; j < dim2; ++j) {
            dptr[i][j] = new double[dim3];
        }
    }
    #pragma omp critical
    {
        
    }
    
}

//functions in verlet.py
void MD::verletNextR(double **& nR, double ** r_t, double ** v_t, double ** acc_t, double h){
    #pragma omp parallel for
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            nR[i_atom][i_dim] = r_t[i_atom][i_dim] + h * v_t[i_atom][i_dim] + 0.5 * h * h * acc_t[i_atom][i_dim];
        }
    }
    #pragma omp critical
    {
        
    }
}

void MD :: verletNextV(double **& nV, double** v_t, double** acc_t, double** acc_t_plus_h, double h){
    #pragma omp parallel for
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            nV[i_atom][i_dim] = v_t[i_atom][i_dim] + 0.5 * h * (acc_t[i_atom][i_dim] + acc_t_plus_h[i_atom][i_dim]) ;
        }
    }
    #pragma omp critical
    {
        
    }
}


//functions in properties.py
void MD::my_distance(double *& drij){
    MD::my_distance_ = 0;
    for (int i = 0; i < dim; ++i) {
        MD::my_distance_ += drij[i]*drij[i];
    }
    MD::my_distance_ = pow(MD::my_distance_,0.5);
}

void MD::my_disp_in_box(double *& drij){
    for (int i = 0; i < dim; ++i) {
        drij[i] -= L * round(drij[i]/L);
    }
    
    
}

void MD::my_pos_in_box(double **& pos){
    #pragma omp parallel for
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            pos[i_atom][i_dim] -= L * round(pos[i_atom][i_dim]/L);
        }
    }
    #pragma omp critical
    {
        
    }

}

void MD::my_kinetic_energy(double **& vel){
    my_kinetic_energy_ = 0;
    #pragma omp parallel for
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
             my_kinetic_energy_ += vel[i_atom][i_dim] * vel[i_atom][i_dim];
        }
    }
    my_kinetic_energy_ *= 0.5 * M ;
    #pragma omp critical
    {
        
    }


}

void MD::get_displacement_table(double **& R){
    double * disp = new double[dim];
    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1; i_atom2 < N; ++i_atom2) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                disp[i_dim] = R[i_atom1][i_dim] - R[i_atom2][i_dim];
            }
            my_disp_in_box(disp);
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                my_displacement_table_[i_atom1][i_atom2][i_dim] =  disp[i_dim];
                my_displacement_table_[i_atom2][i_atom1][i_dim] =  disp[i_dim];
            }
             
        }
    }
    #pragma omp critical
    {
        
    }

    
}

void MD::get_distance_table(double ***& disp_table){
    double * disp = new double[dim];
    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1; i_atom2 < N; ++i_atom2) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                disp[i_dim] = disp_table[i_atom1][i_atom2][i_dim];
            }
            my_distance(disp);
            my_distance_table_[i_atom1][i_atom2] = my_distance_;
            my_distance_table_[i_atom2][i_atom1] = my_distance_;
             
        }
    }
    #pragma omp critical
    {
        
    }

}

void MD::my_potential_energy(double **& dist_table){
    double vshift = 4 * pow(rc, -6) * (pow(rc,-6)-1);
    my_potential_energy_ = 0;
    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1 + 1; i_atom2 < N; ++i_atom2) {
            double r = dist_table[i_atom1][i_atom2];
            if (r <= rc){
                my_potential_energy_ += 4 * pow(r, -6) * (pow(r,-6)-1) - vshift;
                
            }
            
        }
    }
    #pragma omp critical
    {
        
    }


}

void MD::my_force_on(int i_tagged, double **& pos ){
    #pragma omp parallel for
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
        my_force_on_[i_dim] = 0; // initialize the force to zero everytime
    }
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        if (pos[i_atom][0] == pos [i_tagged][0] && pos[i_atom][1] == pos [i_tagged][1] && pos[i_atom][2] == pos [i_tagged][2]) {
            continue;
        }
        double * rij = new double[dim];
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            rij[i_dim] = pos[i_atom][i_dim] - pos [i_tagged][i_dim];
            
        }
        my_distance(rij);
        if (my_distance_ <= rc) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                my_force_on_[i_dim] += 24 * pow(my_distance_,-8) * (2 * pow(my_distance_,-6) - 1) * rij[i_dim] ;
            }
        }
    }
    #pragma omp critical
    {
        
    }
}

void MD::my_temperature(double & _my_kinetic_energy_){
    my_temperature_ = _my_kinetic_energy_ / (3 * N / 2.0);
}

void MD::my_pressure(double & T, double **& R, double **& Forces){
    double virial = 0;
    double V = pow(L,3);
    for (int i_atom = 0; i_atom < N; ++ i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            virial += R[i_atom][i_dim] * Forces[i_atom][i_dim];
        }
    }
    my_pressure_ = (N * T + (virial /3)) / V;
}
