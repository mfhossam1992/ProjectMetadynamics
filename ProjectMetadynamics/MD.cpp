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
MD::MD(Init*& init_, bool anderson_, double Ta_, double eta_,bool mtd_, double rc_, double metarc_, double h_, string outputfileName, int steps_, string trajFileName_):
    N(init_->getN()),
    dim(init_->getDim()),
    L(init_->getL()),
    T0(init_->getT()),
    anderson(anderson_),
    Ta(Ta_),
    eta(eta_),
    mtd(mtd_),
    rc(rc_),
    metarc(metarc_),
    h(h_),
    output_fileName(outputfileName),
    steps(steps_),
    traj_filename(trajFileName_)

{
    alloc_mem(N, dim);
    alloc_mem_rv(N, dim);
    //double ** Rptr = init_->getPosition();
    //double ** Vptr = init_->getVelocity();
    //for (int i_atom = 0; i_atom < N; ++i_atom) {
    //    for (int i_dim; i_dim < dim; ++i_dim) {
    //        R[i_atom][i_dim] = Rptr[i_atom][i_dim];
    //        V[i_atom][i_dim] = Vptr[i_atom][i_dim];
    //    }
    //}
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
//#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        R[i] = new double[dim_]; // 3-dim
        V[i] = new double [dim_]; // 3-dim
    }
//#pragma omp critical
    {
        
    }
}
//memory allocation for new position and new velocity
void MD::alloc_mem_rv(int N_, int dim_=3){
    nR = new double *[N_];
    nV = new double *[N_];
//#pragma omp parallel for
    for (int i =0; i<N_; ++i) {
        nR[i] = new double[dim_]; // 3-dim
        nV[i] = new double [dim_]; // 3-dim
    }
//#pragma omp critical
    {
        
    }
}

//memory allocation for passed double pointer
void MD::alloc_mem2(double ** & dptr, int dim1, int dim2){
    dptr = new double *[dim1];
//#pragma omp parallel for
    for (int i =0; i<dim1; ++i) {
        dptr[i] = new double[dim2];
    }
//    #pragma omp critical
    {
        
    }
    
}
//memory allocation for passed triple pointer
void MD::alloc_mem3(double *** & dptr, int dim1, int dim2, int dim3){
    dptr = new double **[dim1];
//#pragma omp parallel for
    for (int i =0; i<dim1; ++i) {
        dptr[i] = new double * [dim2];
        for (int j = 0; j < dim2; ++j) {
            dptr[i][j] = new double[dim3];
        }
    }
//    #pragma omp critical
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
//    #pragma omp parallel for
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
//    #pragma omp critical
    {
        
    }

    
}

void MD::get_distance_table(double ***& disp_table){
    double * disp = new double[dim];
//    #pragma omp parallel for
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
//    #pragma omp critical
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
            rij[i_dim] = pos[i_tagged][i_dim] - pos [i_atom][i_dim];
            my_disp_in_box(rij);
            
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

// output.py
void MD::output(string fileName, string line){
    cout << line << endl;
    ofstream output_file;
    output_file.open(fileName, ios_base::app);
    output_file << line << endl;
    output_file.close();
    
}

void MD::simulate(){
    MD::output(output_fileName, "steps, temperature, pressure, energy, Q6\n");
    alloc_mem2(F, N, dim);
    alloc_mem2(A, N, dim);
    alloc_mem2(nF, N, dim);
    alloc_mem2(nA, N, dim);
    
    for (int i_t = 0; i_t < steps; ++i_t) {
        //FOR DEBUG
        //cout <<"\n\n\n\n\n\n\n THIS IS THE "<<to_string(i_t)<<"TH ITERATION\n\n\n\n\n\n\n\n";
        //END FOR DEBUG
        // Anderson Thermostat
        if (anderson == true) {
            double sigma = pow((Ta / M), 0.5);
            double mean = 0;
            for (int i_atom = 0; i_atom < N; ++i_atom) {
                if ((rand()/double(RAND_MAX)) < eta * h) {
                    for (int i_dim = 0; i_dim < dim; ++i_dim) {
                        std::random_device rd{};
                        std::mt19937 gen{rd()};
                        std::normal_distribution<> dist{mean,sigma};
                        V[i_atom][i_dim] = dist(gen);
                    }

                }
            }
        }

        // *********** Propagation ************** //
        // Calculate Forces
        #pragma omp parallel for
        for (int i_atom = 0; i_atom < N; ++i_atom) {
            my_force_on(i_atom, R);
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                F[i_atom][i_dim] = my_force_on_[i_dim];
                A[i_atom][i_dim] = (1/M) * F[i_atom][i_dim];
            }
        }
        #pragma omp critical
        {
            
        }
        // Calculate new Positions
        verletNextR(nR, R, V, A, h);
//        //DEBUG
//        cout << "\n Unwrapped NEW POSITION\n" ;
//        for (int i_atom = 0; i_atom < N; ++i_atom) {
//            cout << nR[i_atom][0] << " , " << nR[i_atom][1] << " , " << nR[i_atom][2] << endl;
//
//        }
        my_pos_in_box(nR);
//        cout << "\n NEW POSITION\n" ;
//        for (int i_atom = 0; i_atom < N; ++i_atom) {
//            cout << nR[i_atom][0] << " , " << nR[i_atom][1] << " , " << nR[i_atom][2] << endl;
//
//        }
//
//        //END_DEBUG

        // Calculate Forces with new Positions
        #pragma omp parallel for
        for (int i_atom = 0; i_atom < N; ++i_atom) {
            my_force_on(i_atom, nR);
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                nF[i_atom][i_dim] = my_force_on_[i_dim];
                nA[i_atom][i_dim] = (1/M) * nF[i_atom][i_dim];
            }
        }
        #pragma omp critical
    {
            
        }
        // Get Displacement Table
        get_displacement_table(nR);
        get_distance_table(my_displacement_table_);
        // MetaDynamics Bias Forces (Still want to implement)
        
        // Calculate New velocities
        verletNextV(nV, V, A, nA, h);
        // update positions
#pragma omp parallel for
        for (int i_atom = 0; i_atom < N; ++i_atom) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                R[i_atom][i_dim] = nR [i_atom][i_dim];
                V[i_atom][i_dim] = nV [i_atom][i_dim];
            }
        }
        #pragma omp critical
        {
            
        }
        // Measuring physical Quantities
        my_kinetic_energy(V);
        my_temperature(my_kinetic_energy_);
        my_potential_energy(my_distance_table_);
        double E_tot = my_kinetic_energy_ + my_potential_energy_;
        my_pressure(my_temperature_, R, nF);
                // left is meta_Q6 method -- calculate_Q6()
        
        // output results
        if (i_t  % 50 == 0) {
            int step = i_t + 1;
            //cout.precision(15);
            //cout<<my_temperature_;
            #pragma omp critical
            {
                
            string output_line = to_string(step) + "  " + to_string(my_temperature_) + "  " + to_string(my_pressure_) + "  " + to_string(E_tot)+"  "; // remaining + Q6
            output(output_fileName, output_line);
            write_xyz(traj_filename, R);
        }
            
        }
    }
    
}

void MD::write_xyz(string traj_file_name, double ** & pos){
    ofstream output_file;
    output_file.open(traj_file_name, ios_base::app);
    output_file << to_string(N) << endl;
    output_file << "Atom Positions" << endl;
    for (int i_atom = 0; i_atom < N; ++i_atom) {
           output_file << "Ar\t" << to_string(pos[i_atom][0])<<"  "<< to_string(pos[i_atom][1])<<"  "<<to_string(pos[i_atom][2])<<"  "<<endl;
            }
    output_file.close();
    
}
