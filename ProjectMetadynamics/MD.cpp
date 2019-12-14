//
//  MD.cpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/8/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#include "MD.hpp"

using namespace std;

// Constructor1_wo_mtd_bias
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
    meta_rc(metarc_),
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
    if (mtd_ == true) {
        alloc_mem2(ds_dr, N, dim);
        alloc_mem2(meta_nR, N, dim);
        alloc_mem3(meta_n_my_displacement_table_, N, N, dim);
        alloc_mem2(meta_n_my_distance_table_, N, N);

    }
    
}
// Constructor2_w_mtd_bias
MD::MD(Init*& init_, bool anderson_, double Ta_, double eta_,bool mtd_, double meta_w_, double meta_sig_, double meta_sig_2_,int max_n_gauss_, int meta_tau_, double rc_, double metarc_, double h_, string outputfileName, int steps_, string trajFileName_):
    N(init_->getN()),
    dim(init_->getDim()),
    L(init_->getL()),
    T0(init_->getT()),
    anderson(anderson_),
    Ta(Ta_),
    eta(eta_),
    mtd(mtd_),
    meta_w(meta_w_),
    meta_sig(meta_sig_),
    meta_sig_2(meta_sig_2_),
    max_n_gauss(max_n_gauss_),
    meta_tau(meta_tau_),
    rc(rc_),
    meta_rc(metarc_),
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
    if (mtd_ == true) {
        alloc_mem2(ds_dr, N, dim);
        alloc_mem2(meta_nR, N, dim);
        alloc_mem3(meta_n_my_displacement_table_, N, N, dim);
        alloc_mem2(meta_n_my_distance_table_, N, N);
        n_gauss = 0;

    }
    
}

// Constructor1'_wo_mtd_bias_post_equil
MD::MD(MD*& init_, bool anderson_, double Ta_, double eta_,bool mtd_, double rc_, double metarc_, double h_, string outputfileName, int steps_, string trajFileName_):
    N(init_->getN()),
    dim(init_->getDim()),
    L(init_->getL()),
    T0(init_->getT()),
    anderson(anderson_),
    Ta(Ta_),
    eta(eta_),
    mtd(mtd_),
    rc(rc_),
    meta_rc(metarc_),
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
    if (mtd_ == true) {
        alloc_mem2(ds_dr, N, dim);
        alloc_mem2(meta_nR, N, dim);
        alloc_mem3(meta_n_my_displacement_table_, N, N, dim);
        alloc_mem2(meta_n_my_distance_table_, N, N);

    }
    
}
// Constructor2'_w_mtd_bias_post_equil
MD::MD(MD*& init_, bool anderson_, double Ta_, double eta_,bool mtd_, double meta_w_, double meta_sig_, double meta_sig_2_,int max_n_gauss_, int meta_tau_, double rc_, double metarc_, double h_, string outputfileName, int steps_, string trajFileName_):
    N(init_->getN()),
    dim(init_->getDim()),
    L(init_->getL()),
    T0(init_->getT()),
    anderson(anderson_),
    Ta(Ta_),
    eta(eta_),
    mtd(mtd_),
    meta_w(meta_w_),
    meta_sig(meta_sig_),
    meta_sig_2(meta_sig_2_),
    max_n_gauss(max_n_gauss_),
    meta_tau(meta_tau_),
    rc(rc_),
    meta_rc(metarc_),
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
    if (mtd_ == true) {
        alloc_mem2(ds_dr, N, dim);
        alloc_mem2(meta_nR, N, dim);
        alloc_mem3(meta_n_my_displacement_table_, N, N, dim);
        alloc_mem2(meta_n_my_distance_table_, N, N);
        n_gauss = 0;

    }
    
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

void MD::get_displacement_table(double *** & disp_table, double **& R){
    double * disp = new double[dim];
//    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1; i_atom2 < N; ++i_atom2) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                disp[i_dim] = R[i_atom1][i_dim] - R[i_atom2][i_dim];
            }
            my_disp_in_box(disp);
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                disp_table[i_atom1][i_atom2][i_dim] =  disp[i_dim];
                disp_table[i_atom2][i_atom1][i_dim] =  disp[i_dim];
            }
             
        }
    }
//    #pragma omp critical
    {
        
    }
    delete [] disp;
    
}

void MD::get_distance_table(double ** & dist_table, double ***& disp_table){
    double * disp = new double[dim];
//    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1; i_atom2 < N; ++i_atom2) {
            for (int i_dim = 0; i_dim < dim; ++i_dim) {
                disp[i_dim] = disp_table[i_atom1][i_atom2][i_dim];
            }
            my_distance(disp);
            dist_table[i_atom1][i_atom2] = my_distance_;
            dist_table[i_atom2][i_atom1] = my_distance_;
             
        }
    }
//    #pragma omp critical
    {
        
    }
    delete [] disp;

}

void MD::my_potential_energy(double **& dist_table){
    double vshift = 4 * pow(rc, -6) * (pow(rc,-6)-1);
    my_potential_energy_ = 0;
//    #pragma omp parallel for
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1 + 1; i_atom2 < N; ++i_atom2) {
            double r = dist_table[i_atom1][i_atom2];
            if (r <= rc){
                my_potential_energy_ += 4 * pow(r, -6) * (pow(r,-6)-1) - vshift;
                
            }
            
        }
    }
//    #pragma omp critical
//    {
//
//    }


}

void MD::my_force_on(int i_tagged, double **& pos ){
    double * rij = new double[dim];
   #pragma omp parallel for
    for (int i_dim = 0; i_dim < dim; ++i_dim) {
        my_force_on_[i_dim] = 0; // initialize the force to zero everytime
    }
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        if (pos[i_atom][0] == pos [i_tagged][0] && pos[i_atom][1] == pos [i_tagged][1] && pos[i_atom][2] == pos [i_tagged][2]) {
            continue;
        }
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
    delete [] rij;
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
    if (meta_sig_2 == 0) {
        MD::output(output_fileName, "steps, temperature, pressure, energy, Q6\n");

    } else {
        MD::output(output_fileName, "steps, temperature, pressure, energy, Q6, PE\n");

    }
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
            //modification for temperature blowup problem
            double my_eta_ = eta;
            if (my_temperature_ > 2.5) {
                my_eta_ = (my_temperature_ / 2 )* eta; // higher thermostat coupling !!
            }
            for (int i_atom = 0; i_atom < N; ++i_atom) {
                if ((rand()/double(RAND_MAX)) < my_eta_ * h) {
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
        get_displacement_table(my_displacement_table_,nR);
        get_distance_table(my_distance_table_,my_displacement_table_);
        my_potential_energy(my_distance_table_); // modified its position to be used in metaD_bias

        // MetaDynamics Bias Forces (Still want to implement)
        if (mtd == true) {
            if (meta_sig_2 != 0) {
                meta_2(nF, i_t, nR);
            }
            else meta(nF, i_t, nR);
                  
            for (int i_atom = 0; i_atom < N; ++i_atom) {
                for (int i_dim = 0; i_dim < dim; ++i_dim) {
                    nA[i_atom][i_dim] = (1/M) * nF[i_atom][i_dim];
                }
            }
        }
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
        //my_potential_energy(my_distance_table_);
        double E_tot = my_kinetic_energy_ + my_potential_energy_; // my_potential_energy_ is calculated before MTD BIAS because we use it in bias as well
        my_pressure(my_temperature_, R, nF);
                // left is meta_Q6 method -- calculate_Q6()
        
        // output results
        //if (i_t  % 50 == 0) {
            //int step = i_t + 1;
            //cout.precision(15);
            //cout<<my_temperature_;
         //   #pragma omp critical
         //   {
                if (mtd == false) {
                    calculate_Q6(Q_6,my_distance_table_, my_displacement_table_);
                    string output_line = to_string(i_t) + "  " + to_string(my_temperature_) + "  " + to_string(my_pressure_) + "  " + to_string(E_tot)+"  "+ to_string(Q_6);
                    output(output_fileName, output_line);



                }
                else{
                    if (meta_sig_2 != 0) {
                        string output_line = to_string(i_t) + "  " + to_string(my_temperature_) + "  " + to_string(my_pressure_) + "  " + to_string(E_tot)+"  "+ to_string(meta_Q6)+"  "+ to_string(my_potential_energy_);
                        output(output_fileName, output_line);
                    }
                    else {
                    string output_line = to_string(i_t) + "  " + to_string(my_temperature_) + "  " + to_string(my_pressure_) + "  " + to_string(E_tot)+"  "+ to_string(meta_Q6);
                    output(output_fileName, output_line);
                    }

                }
        if (i_t % 100 == 0) {
                        write_xyz(traj_filename, R);

        }
       // }
            
       // }
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


void MD::calculate_Qlm(int l, int m, double ** & dist_table, double *** & displacement_table){
    Q_lm_r_ = 0; // initialize to zero
    Q_lm_i_ = 0; // initialize to zero
    double n_bonds = 0;
    double * rhat = new double[3];
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_atom2 = i_atom1 + 1; i_atom2 < N; ++ i_atom2) {
            double r = dist_table[i_atom1][i_atom2];
            if (r <= meta_rc) {
                n_bonds += 1;
                for (int i_dim = 0; i_dim < dim; ++i_dim) {
                    rhat[i_dim] = displacement_table[i_atom1][i_atom2][i_dim] / r;
                }
                double theta = acos(rhat[2]);
                double phi = atan2(rhat[1], rhat[0]);
                Q_lm_r_ += boost::math::spherical_harmonic_r(l, m, theta, phi);
                Q_lm_i_ += boost::math::spherical_harmonic_i(l, m, theta, phi);

            }

        }
    }
    if (n_bonds != 0) {
        Q_lm_r_ /= n_bonds;
        Q_lm_i_ /= n_bonds;
    }
    delete [] rhat;
    
}

void MD::calculate_Ql(int l, double ** & dist_table, double *** & displacement_table){
    Q_l = 0; // initialize to zero
    calculate_Qlm(l,0,dist_table,displacement_table); // calculate m = 0 term first and append it
    Q_l = pow(Q_lm_r_, 2) + pow(Q_lm_i_, 2);
    for (int i_m = 1; i_m < (l+1); ++i_m) { // range of i_m is from 1 to l instead from -l to l because of symmetry around i_m = 0
        calculate_Qlm(l,i_m,dist_table,displacement_table);
        Q_l += 2 * (pow(Q_lm_r_, 2) + pow(Q_lm_i_, 2)); // because the squares are the same from i_m -l range as the +l range, excluding the m = zero term, then I calculate one side and multiply by 2 instead of going from -l to l in the loop
        
    }
    Q_l = sqrt((4 * M_PI * Q_l)/(2 * l + 1));
}

void MD::calculate_Q6(double & Q6, double ** & dist_table_, double *** & displacement_table_){
    calculate_Ql(6,dist_table_, displacement_table_);
    Q6 = Q_l;
}

void MD::calculate_ds_dr(double **  & pos){
    double d = 0.001;
    //make copies of current pos, disp_table, distance_table and store them
    for (int i_atom = 0; i_atom < N; ++i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            meta_nR[i_atom][i_dim] = pos[i_atom][i_dim];
        }
    }
    get_displacement_table(meta_n_my_displacement_table_, meta_nR); // makes a copy of current displacement table by actually re-calculate it since it is of comparable complexity
    get_distance_table(meta_n_my_distance_table_, meta_n_my_displacement_table_);//makes a copy of current displacement table by actually re-calculate it since it is of comparable complexity
    calculate_Q6(meta_Q6, my_distance_table_, my_displacement_table_); // calculates the current un-shifted value of Q_6
    // AS for my_potential_energy_ it is already calculated in Simulate function just before calling the meta function to bias forces
    // END_OF_MAKING_COPIES_AND_CALCULATIONS_Prior_TO_SHIFTING
    
    //START SHIFTING
    double * disp_ = new double[dim];
    for (int i_atom1 = 0; i_atom1 < N; ++i_atom1) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            meta_nR[i_atom1][i_dim] += d;
            for (int i_atom2 = 0; i_atom2 < N; ++i_atom2) {
                for (int i_dim_ = 0; i_dim_ < dim; ++ i_dim_) {
                    disp_[i_dim_] = meta_nR[i_atom1][i_dim_] - meta_nR[i_atom2][i_dim_];
                }
                my_disp_in_box(disp_);
                for (int i_dim_ = 0; i_dim_ < dim; ++ i_dim_) {
                    meta_n_my_displacement_table_[i_atom1][i_atom2][i_dim_] = disp_[i_dim_];
                    meta_n_my_displacement_table_[i_atom2][i_atom1][i_dim_] = disp_[i_dim_];
                }
                my_distance(disp_);
                meta_n_my_distance_table_[i_atom1][i_atom2] = my_distance_;
                meta_n_my_distance_table_[i_atom2][i_atom1] = my_distance_;
                
            }
            calculate_Q6(meta_n_Q6, meta_n_my_distance_table_, meta_n_my_displacement_table_);
            ds_dr[i_atom1][i_dim] = (meta_n_Q6 - meta_Q6) / d;
            meta_nR[i_atom1][i_dim] -= d;
        }
    }
    delete [] disp_;
    
}


void MD::meta(double ** & nF_, int i_t_, double ** & pos){
   // number of gaussians by ref(to be incremented continuously, vector of updated gaussian center positions,
    // calculate the derivative of s w.r.t atom positions
    calculate_ds_dr(pos);
    // every tau step, save the value of s = meta_Q6
    if (i_t_ % meta_tau == 0) {
        n_gauss += 1;
        if (n_gauss < max_n_gauss) {
            S.push_back(meta_Q6);
            string line = to_string(i_t_) + "  " + to_string(meta_Q6);
            output("Q6_Gauss_Center.txt", line);
            
        } else {
            cout<<"MAX_NUMBER_OF_GAUSS_EXCEEDED"<<endl;
            exit(EXIT_SUCCESS);
        }
    }
    //calculate the derivative of the history-dependent potential w.r.t s=meta_Q6
    double dV_ds = 0;
    for (auto & s_tau : S) {
        //DEBUG
        //cout << s_tau<<" ; "<< meta_Q6 << endl;
        // END_DEBUG
        double gauss = meta_w * exp(- pow((meta_Q6 - s_tau), 2) / (2 * pow(meta_sig, 2)));
        dV_ds += (-1) * gauss * (meta_Q6 - s_tau) / pow(meta_sig, 2);
    }
    // Bias_the_force
    for (int i_atom = 0; i_atom < N; ++ i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            nF_[i_atom][i_dim] += (-1) * dV_ds * ds_dr[i_atom][i_dim];
        }
    }
}

void MD::meta_2(double ** & nF_, int i_t_, double ** & pos){
    // calculate the derivative of s w.r.t atom positions -- this also calculates inside it the current value of meta_Q6
    calculate_ds_dr(pos);
    // every tau step, save the value of s = meta_Q6 and s_2 = my_potential_energy_
    if (i_t_ % meta_tau == 0) {
        n_gauss += 1;
        if (n_gauss < max_n_gauss) {
            S.push_back(meta_Q6);
            S_2.push_back(my_potential_energy_);
            string line = to_string(i_t_) + "  " + to_string(meta_Q6)+ "  " + to_string(my_potential_energy_);
            output("Q6_PE_Gauss_Centers.txt", line);
            
        } else {
            cout<<"MAX_NUMBER_OF_GAUSS_EXCEEDED"<<endl;
            exit(EXIT_SUCCESS);
        }
    }
    // Calculate the bias force and add it for each tau step
    size_t N_tau = S.size(); // or N_tau = S_2.size()
    for (int i_atom = 0; i_atom < N; ++ i_atom) {
        for (int i_dim = 0; i_dim < dim; ++i_dim) {
            double f_tau = 0; // the biasing force as per each atom per each dimension
            for (size_t i_tau = 0; i_tau < N_tau; ++ i_tau) {
                double s_tau = S[i_tau];
                double s_2_tau = S_2[i_tau];
                double gauss = meta_w * exp(- pow((meta_Q6 - s_tau), 2) / (2 * pow(meta_sig, 2))) * exp(- pow((my_potential_energy_ - s_2_tau), 2) / (2 * pow(meta_sig_2, 2)));
                f_tau += gauss * ((((meta_Q6 - s_tau)/pow(meta_sig, 2)) * ds_dr[i_atom][i_dim]) + (((my_potential_energy_ - s_2_tau)/pow(meta_sig_2, 2)) * (-1) * (nF_[i_atom][i_dim])));
            }
            nF_[i_atom][i_dim] += f_tau;
        }
    }


}
//geters
double ** MD::getPosition(){
    return MD::R;
}

double ** MD::getVelocity(){
    return MD::V;
}
int MD::getN(){
    return MD::N;
}
double MD::getT(){
    return MD::my_temperature_;
}
double MD::getL(){
    return MD::L;
}
int MD::getDim(){
    return MD::dim;
}

//*end geters*//
