//
//  MD.hpp
//  ProjectMetadynamics
//
//  Created by Hossam Farag on 12/8/19.
//  Copyright © 2019 Hossam Farag. All rights reserved.
//

#ifndef MD_hpp
#define MD_hpp

#include <stdio.h>
#include "Init.hpp"
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <random>
#include <cstring>

//#ifdef BOOST
#include <boost/math/special_functions/spherical_harmonic.hpp>
//#endif
using namespace std;

class MD {
    //private variables
    double ** R;
    double ** V;
    double ** nR;
    double ** nV;
    double ** F;
    double ** A;
    double ** nF;
    double ** nA;
    vector<double> S; // updated gaussian center positions for Q6
    vector<double> S_2; // updated gaussian center positions for my_potential_energy_
    //int n_gauss = 0; // number of added gaussians (Re-Decalred beneath with bunch of other MTD variables)
    double Ta; // Anderson theromstat temperature
    double eta; // Anderson thermostat parameter
    string output_fileName;
    string traj_filename;
    string traj_v_filename;
    int steps;
    //same variables as in Init Class
    int N ;
    int dim = 3; // default
    double L = 4; // default = 4
    double T0 = 1; // initial Temperature (from init->getT())
    double M = 48; // default for Argon 
    // gain from constructor
    bool anderson;//do thermostat or not
    bool mtd; // perform mtd or not
    double rc;//LJ cutoff
    double meta_rc;// n.n. cutoff for Q6
    double h;//timestep
    //variables corresponding to the return of each function which became void
    double my_distance_ = 0;//variable updated via my_distance() method
    double my_kinetic_energy_ = 0;//updated via my_kinetic_energy() method
    double *** my_displacement_table_;//(natom,natom,ndim) displacement table updated by get_displacement_table() method when passedto it
    double ** my_distance_table_;//(natom,natom) distance table updated by get_distance_table() method hen passed to it
    double my_potential_energy_;// updated by my_potential_energy() method
    double * my_force_on_; // updated by my_force_on() method
    double my_temperature_; //updated by my_temperature() method
    double my_pressure_; //updated by my_pressure() method
    double Q_lm_r_; // updated by calculate_Qlm() method
    double Q_lm_i_; // updated by calculate_Qlm () method
    double Q_l; // updated by calculate_Ql() method
    double Q_6; // updated by calculate_Q6() method
    //variables corresponding to metadynamics bias
    double ** ds_dr; // (natom,ndim) the derivative of s with respect to atom positions
    double ** meta_nR; // copy of current R to shift it and do finite difference
    double *** meta_n_my_displacement_table_; // displacement table corresponding to the shifted positions
    double ** meta_n_my_distance_table_; // distance table corresponding to the shifted positions
    double meta_Q6; // Q6 corresponding to positions before shifting
    double meta_n_Q6; //Q6 corresponding to shifted positions
    int n_gauss; // keeping track of the number of added gaussians and terminate the program when this number exceeds max_n_gauss
    // MTD_Gaussian_Parameters
    double meta_w;
    double meta_sig;
    double meta_sig_2 = 0;
    int max_n_gauss;
    int meta_tau;
    //private methods
    // memory allocation
    void alloc_mem(int,int);
    void alloc_mem_rv(int, int);
    void alloc_mem2(double ** &,int, int);
    void alloc_mem3(double *** &, int, int, int);
    //functions in verlet.py
    void verletNextR (double ** &, double **, double **, double **, double );//(natom,ndim) next position array, (natom,ndim) position array, (natom,ndim) velocity array, (natom,ndim) acceleration array,time step h
    void verletNextV (double **&, double**,double**,double**,double);//(natom,ndim) next velocity array, (natom,ndim) velocity array,(natom,ndim) acceleration array, (natom,ndim) next time step acceleration array,time step h
    
    //functions in properties.py
    void my_distance(double *&); // displacement vector of length 3 assuming it is already accounted for PBC
    void my_disp_in_box(double *& ); // length-3 displacement vector
    void my_pos_in_box (double **&); // (natom,ndim) position array
    void my_kinetic_energy (double **& ); //(natom,ndim) velocity array
    void get_displacement_table (double *** &, double **&); // the (natom,natom,ndim) varibale the function will update, (natom,ndim) position array
    void get_distance_table(double ** &, double***& );//the (natom,natom) varibale the function will update, (natom,natom,ndim) displacement table
    void my_potential_energy(double **&); // (natom,natom) distance table
    void my_force_on(int, double **& ); // particle index,(natom,ndim) position array
    void my_temperature(double&); //kinetic energy
    void my_pressure( double&, double **&,double **&); // Temperature,(natom,ndim) position array,(natom,ndim) forces array
    
    //functions in metadynamics.py
    void calculate_Qlm(int, int, double ** &, double *** &); //order l, order -l<m<l, distance table, displacement table
    void calculate_Ql(int, double ** &, double *** &); // order l, meta_rc, distance table, displacement table
    void calculate_Q6(double &, double ** &, double *** &); // Q_6 or meta_n_Q6 or meta_Q6, distance table, displacement table, runs calculate_Ql with order 6 and assign the value of Q_l to Q_6 or meta_n_Q6 or meta_Q6
    void calculate_ds_dr(double ** &); //(natom,ndim) position array, corresponding Q6 to position array
    void meta(double ** &, int, double ** &);//(natom,ndim) Force varibale the function will update, step number,(natom,ndim) position array
    void meta_2(double ** &, int, double ** &); //(natom,ndim) Force varibale the function will update, step number,(natom,ndim) position array
    
    //functions in output.py
    void output(string, string); // filename, thermo_output
    void write_xyz(string, double ** &); //filename,(natom,ndim) position array
public:
    MD(Init*&, bool anderson, double Ta_, double eta_,bool mtd, double rc_, double meta_rc_, double h_, string output_fileName_, int steps,string trajFileName); // pointer to Init Object, anderson, thermostat desired temperature, thermostat parameter,mtd, lj rc, mtd nn cutoff, timestep size, outputfileName, total number of MD steps, xyz trajectory file name
//    MD(Init*&, bool anderson, double Ta_, double eta_,bool mtd, double meta_w_, double meta_sig_, int max_n_gauss_, int meta_tau_, double rc_, double meta_rc_, double h_, string output_fileName_, int steps,string trajFileName); // pointer to Init Object, anderson, thermostat desired temperature, thermostat parameter,mtd, mtd_gauss_height, mtd_gauss_width , maximum_num_gauss, gauss_deposition_frequency,lj rc, mtd nn cutoff, timestep size, outputfileName, total number of MD steps, xyz trajectory file name
    MD(Init*&, bool anderson, double Ta_, double eta_,bool mtd, double meta_w_, double meta_sig_, double meta_sig_2,int max_n_gauss_, int meta_tau_, double rc_, double meta_rc_, double h_, string output_fileName_, int steps,string trajFileName); // pointer to Init Object, anderson, thermostat desired temperature, thermostat parameter,mtd, mtd_gauss_height, mtd_gauss_width , maximum_num_gauss, gauss_deposition_frequency,lj rc, mtd nn cutoff, timestep size, outputfileName, total number of MD steps, xyz trajectory file name
    MD(MD*&, bool anderson, double Ta_, double eta_,bool mtd, double rc_, double meta_rc_, double h_, string output_fileName_, int steps,string trajFileName); // pointer to Init Object, anderson, thermostat desired temperature, thermostat parameter,mtd, lj rc, mtd nn cutoff, timestep size, outputfileName, total number of MD steps, xyz trajectory file name
    MD(MD*&, bool anderson, double Ta_, double eta_,bool mtd, double meta_w_, double meta_sig_, double meta_sig_2,int max_n_gauss_, int meta_tau_, double rc_, double meta_rc_, double h_, string output_fileName_, int steps,string trajFileName); // pointer to Init Object, anderson, thermostat desired temperature, thermostat parameter,mtd, mtd_gauss_height, mtd_gauss_width , maximum_num_gauss, gauss_deposition_frequency,lj rc, mtd nn cutoff, timestep size, outputfileName, total number of MD steps, xyz trajectory file name

    

    virtual ~MD();
    void simulate(); // main function to perform simulation
    void equilibrate(); // main function to perform simulation
    //double * simulate_mtd();
    //getters
    double ** getPosition();
    double ** getVelocity();
    int getN();
    double getT();
    double getL();
    int getDim();
    //Setters
    void setTa(double & Ta_); //sets the thermostat Temperature
   

    
    
};
#endif /* MD_hpp */
