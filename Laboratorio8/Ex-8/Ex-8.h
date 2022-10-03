#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.h"
#include <vector>
int seed[4];
Random rnd;

int SA;

// simulation
int nstep, nblk, stepSA,restart;

//8.1
double span;        //span = x0-xf

// thermodynamical state
double beta,temp,db,temp_o;

//Schrodinger equation
double x_i;
double x0,xf,mu,sigma,d_mu,d_sigma,span_mu,span_sigma,step_dx;
double x_old, x_new;
std :: vector<double> f_schr;                 // array che contiene la definzione della funzione di schrodinger
std :: vector<double> f_schr2;

// Energia
double deltaH;
std :: vector<double> walker;
double blkH_ave,blkH_norm;
double globH_ave,globH_ave2;
double stima_H;
double err_H;
double accepted,attempted;

double x_in,x_fin;      //posizione prima e dopo lo step di montecarlo



//parameters, observables
const int m_props=1000;
int n_props,iu,ic,im,ix,ig;
double nbins;
//double walker[m_props];



//configuration
const int m_spin=50;
double s[m_spin];



//functions
void Input();
void Averages(int);
void Move_MC_metro();
double Eval_H(double);
double Gaussiana_combinaz(double);
double Error(double,double,int);

#endif
