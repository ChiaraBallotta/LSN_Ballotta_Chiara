/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

std :: string esercizio;
// Simulazione NVT - valori delta
double delta_solido = 0.11;
double delta_liquido = 0.2;
double delta_gas = 6;

//radial distribution
const int nbins = 100;
double bin_size;
int gdr[nbins];
double blk_av_gdr[nbins],blk_norm_gdr,glob_av_gdr[nbins],glob_av2_gdr[nbins],err_gdr[nbins];  


//parameters, observables
const int m_props=1000;
int n_props, iv, ik, it, ie, iw;
double vtail, ptail, sd;
double walker[m_props];
   
// averages
double blk_av[m_props], blk_norm, accepted, attempted;
double ave_acc = 0;
double glob_av[m_props], glob_av2[m_props];
double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp,stima_gdr;
double err_pot, err_press, err_kin, err_etot, err_temp;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta,temp,energy,vol,rho,box,rcut;

// simulation
int iNVET, nstep, nblk, restart;
double delta;

//pigreco
const double pi=3.1415927;

//tail corrections
double V_tail,W_tail;

//Ottimizzazione
int step_ottimiz_NVE = 6000; 
int step_ottimiz_NVT = 1000;
int ottimizzazione;

int graf_lim;

//functions
void Input(std :: string);
void Reset(int);
void Accumulate(void);
void Averages(int,std :: string);
void Move(void);
void ConfFinal(std:: string);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);
double tail_corr_W();
double tail_corr_V();
void Ottimizzazione(std::string);

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
