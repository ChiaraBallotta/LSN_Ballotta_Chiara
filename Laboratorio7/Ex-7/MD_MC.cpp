/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <cstring>
#include "MD_MC.h"

using namespace std;

int main()
{ 
  //string Fase = "solid";   // indica su quale stato stiamo svolgendo la simulazione = SOLID, LIQUID,GAS
  //string Fase = "liquid";
  string Fase = "gas";

  Input(Fase); //Inizialization

  //Ottimizzazione
  if(ottimizzazione == 1){
  Ottimizzazione(Fase);
  }

  int nconf = 1;
  for(int iblk=1; iblk <= nblk; iblk++) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; istep++)
    {
      Move();       //implementa l'algoritmo di Verlet
      Measure();    //misura le osservabili durante la simulazione
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk,Fase);   //Print results for current block
  }
  ConfFinal(Fase); //Write final configuration

  cout<<"Acceptance rate media : "<<ave_acc/nblk;

  return 0;
}


void Input(string Fase)
{
  ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  Primes.open("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

//Read input informations
  ReadInput.open("input/input_"+ Fase + ".dat");

  ReadInput >> iNVET;		// decide se fare MD o MC
  if(iNVET == 0) cout<<"Simulazione NVE"<<endl;
  if(iNVET == 1) cout<<"Simulazione NVT"<<endl;

  ReadInput >> restart;		// se c'e 0 riparte dalla configurazione iniziale altrimenti no

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();
// ha finito di settare il generatore
  ReadInput >> ottimizzazione;
  ReadInput >> graf_lim;
  ReadInput >> temp;
  beta = 1.0/temp;	//unita naturali
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;    //legge dove interrompo il potenziale
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  ReadInput >> delta;   //legge il passo di integrazione temporale
  if( iNVET == 1){
    if(Fase == "solid") delta = delta_solido;
    if(Fase == "liquid") delta = delta_liquido;
    if(Fase == "gas") delta = delta_gas;
  }

  ReadInput >> nblk;   

  ReadInput >> nstep;
  bin_size = box/double(nbins * 2);
  cout<<"Dimensione  di bins per g(r) : "<<bin_size<<endl;
  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4; // PRESSIONE
  n_props = 5; //Number of observables

//Read initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart)
  {
    if(iNVET == 0){       //NVE
      ReadConf.open("NVE/output/"+Fase+"/config.out");
      ReadVelocity.open("NVE/output/"+Fase+"/velocity.out");
      for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
    }
    if(iNVET == 1){        //NVT
    ReadConf.open("NVT/output/"+Fase+"/config.out");
    ReadVelocity.open("NVT/output/"+Fase+"/velocity.out");
    for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
    }
  }
  else 
  {
    ReadConf.open("input/config.in");

    cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.Gauss(0.,sqrt(temp));
      vy[i] = rnd.Gauss(0.,sqrt(temp));
      vz[i] = rnd.Gauss(0.,sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    cout << "velocity scale factor: " << fs << endl << endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }
//abbiamo tolto il drift dalle velocita

  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );//periodic boundary conditions
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)   // iNVET == 1 --> NVT:  Si usa montacarlo
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else        // iNVET == 0 --> NVE : algoritmo di Verlet
    {
      xold[i] = Pbc(x[i] - vx[i] * delta);
      yold[i] = Pbc(y[i] - vy[i] * delta);
      zold[i] = Pbc(z[i] - vz[i] * delta);
    }
  }
  
// Calcolo delle TAIL corrections:
   W_tail = tail_corr_W();
   V_tail = tail_corr_V();
//Evaluate properties of the initial configuration
  Measure();

//Print initial values for measured properties
  cout << "Initial potential energy (without tail correction)  = " << walker[iv]/(double)npart << endl;
  cout << "Initial temperature      = " << walker[it] << endl;
  cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << endl;
  cout << "Initial total energy     = " << walker[ie]/(double)npart << endl;
  cout << "Pressure (without tail correction)     = " << walker[iw] << endl;

  return;
}


void Move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.Rannyu()*npart);

    //Old
      energy_old = Boltzmann(x[o],y[o],z[o],o);

    //New
      x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
      y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
      z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

      energy_new = Boltzmann(x[o],y[o],z[o],o);

    //Metropolis test
      p = exp(beta*(energy_old-energy_new));
      if(p >= rnd.Rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i){ //Force acting on particle i
      fx[i] = Force(i,0);
      fy[i] = Force(i,1);
      fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){ //Verlet integration scheme, muoviamo le particelle usando l'algoritmo di Verlet

      xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );  // m = 1
      ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
      znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

      vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];   //contenuto al tempo t - dt
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted = accepted + 1.0;
      attempted = attempted + 1.0;
    }
  }
  return;
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){                       //cut off del potenziale
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Measure() //Properties measurement
{
  double v = 0.0, kin=0.0, w = 0.0;
  double vij,wij;
  double dx, dy, dz, dr;

  for( int i = 0; i < nbins; i++) gdr[i] = 0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
// distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      //Calcolo  g(r)
      for(int b = 0; b < nbins; b++){
        if(dr > (b*bin_size) && dr < ((b+1)*bin_size))  gdr[b] +=2;
      }

      if(dr < rcut)     //condizione che verifica che le particelle siano interagenti
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        v += vij;
	      wij = 48.0*(1.0/pow(dr,12) - 0.5/pow(dr,6));
	      w += wij;
      }
    }          
  }

  for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

// salvo le quantità di interesse e aggiungo le correzzioni di tail a pressione e potenziale
  walker[iv] = 4.0 * v + V_tail; // Potential energy, stima relativa alla particella
  walker[ik] = kin; // Kinetic energy
  walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[iw] = rho*walker[it] + (1.0/ (3* vol))*(w + W_tail);	//pressione

  return;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
       for(int i = 0; i < nbins; i++){
            glob_av_gdr[i] = 0;
            glob_av2_gdr[i] = 0;
            gdr[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
    for(int i = 0; i < nbins; i++){
        blk_av_gdr[i] = 0;
        gdr[i] = 0;
    }
    blk_norm_gdr = 0;
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }

   for(int i = 0; i < nbins; i++){
    blk_av_gdr[i] = blk_av_gdr[i] + gdr[i];
    gdr[i] = 0;
   }
   blk_norm = blk_norm + 1.0;
   blk_norm_gdr = blk_norm_gdr + 1.0;
}


void Averages(int iblk, string Fase) //Print results for current block
{
    
   ofstream Epot, Ekin, Etot, Temp,Press,Gdr,Gdr_err;
   const int wd=18;
    
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
   ave_acc +=  accepted/attempted;
   
   string esercizio;
   if(iNVET == 0) esercizio = "NVE";
   if(iNVET == 1) esercizio = "NVT";

    if(graf_lim== 0){
      Ekin.open(esercizio +"/output/"+Fase+"/output_ekin.dat",ios::app);
      Temp.open(esercizio +"/output/"+Fase+"/output_temp.dat",ios::app);
      Etot.open(esercizio +"/output/"+Fase+"/output_etot.dat",ios::app);  
    }

    Epot.open(esercizio +"/output/"+Fase+ "/output_epot.dat",ios::app);
    Press.open(esercizio +"/output/"+Fase+"/output_press.dat",ios::app);
    Gdr.open(esercizio +"/output/"+Fase+"/rad_distr.dat",ios::app);
    Gdr_err.open(esercizio +"/output/"+Fase+"/rad_distr_err.dat",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    if(graf_lim== 0){
      stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
      glob_av[ik] += stima_kin;
      glob_av2[ik] += stima_kin*stima_kin;
      err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

      stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
      glob_av[ie] += stima_etot;
      glob_av2[ie] += stima_etot*stima_etot;
      err_etot=Error(glob_av[ie],glob_av2[ie],iblk);

      stima_temp = blk_av[it]/blk_norm; //Temperature, non divido per npart siccome è intensiva
      glob_av[it] += stima_temp;
      glob_av2[it] += stima_temp*stima_temp;
      err_temp=Error(glob_av[it],glob_av2[it],iblk);
    }
    stima_pres = blk_av[iw]/blk_norm; //Pressione, non divido per npart siccome è intensiva
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);


  if(graf_lim == 0){
//Kinetic energy
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
  }
//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Pressione
    Press << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_press << endl;

//cout << "----------------------------" << endl << endl;

/////////// 07.3 : calcolo di g(r):
  double r;
  for( int i = 0; i < nbins;i++){
    r = i*bin_size;     // punto d'inizio del bin
    double V_guscio = (4.0/3.0)*rho*M_PI*(pow((r+bin_size),3)-pow(r,3));
    stima_gdr = blk_av_gdr[i]/blk_norm_gdr/V_guscio;
    glob_av_gdr[i] +=stima_gdr;
    glob_av2_gdr[i] +=stima_gdr*stima_gdr;
    err_gdr[i] = Error(glob_av_gdr[i],glob_av2_gdr[i],iblk);
    Gdr<<glob_av_gdr[i]<<"  ";
    Gdr_err<<err_gdr[i]<<"  ";
  }
  Gdr<<endl;
  Gdr_err<<endl;

    Epot.close();
    if(graf_lim == 0){
      Ekin.close();
      Etot.close();
      Temp.close();
    }
    Press.close();
    Gdr.close();
    Gdr_err.close();
}


void ConfFinal(string Fase)
{
  ofstream WriteConf, WriteVelocity, WriteSeed;
  
  string esercizio;
  if(iNVET == 0) esercizio = "NVE";
  if(iNVET == 1) esercizio = "NVT";
  
  cout << "Print final configuration to file config.out" << endl << endl;
  WriteConf.open(esercizio +"/output/"+Fase+"/config.out");
  WriteVelocity.open(esercizio +"/output/"+Fase+"/velocity.out");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
  }
  WriteConf.close();
  WriteVelocity.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

// Calcolo delle tail correction su V
double tail_corr_V(){
  return (8 * M_PI * rho * (1.0/(9*pow(rcut,9))-1.0/(3*pow(rcut,3))) * npart);
}

// Calcolo delle tail correction su P
double tail_corr_W(){
  return (32* M_PI * rho * (1.0/(9*pow(rcut,9))-1.0/(6*pow(rcut,3)))*3*npart);
}


//Esegue l'ottimizzazione
void Ottimizzazione(string Fase){
//Modifico temporaneamente  nstep e nblk
  int nstep_old = nstep;
  int nblk_old = nblk;
  int step = 0;
  if(iNVET == 0)  step = step_ottimiz_NVE;
  if(iNVET == 1)  step = step_ottimiz_NVT;

  for(int i =0; i <step;i++) Move();

//resetto le condizioni di lavoro della simulazione
  nstep = nstep_old;
  nblk = nblk_old;
  cout<<"La fase di ottimizzazione è terminata"<<endl;
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
