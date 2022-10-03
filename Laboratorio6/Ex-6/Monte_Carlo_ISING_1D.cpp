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
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main( int argc, char *argv[])
{  
//tempo a cui valutiamo la simulazione
double temperature = atof(argv[1]);
int camp = atoi(argv[2]);
cout<<"temperatura iniziale:  "<<temperature<<"   "<<argv[2]<<endl;
// Inizializzazione del campionamento scelto
  cout<<"Simulazione del modello di Ising in 1D"<<endl<<"Selezione del metodo di sampling:"<<endl<<endl;
  cout<<"Sampling tramite :  ";
  if(camp == 1) cout<<"METROPOLIS"<<endl;
  if(camp == 2) cout<<"GIBBS"<<endl;
  
  string sim_type;
  if (camp == 1) sim_type = "Metropolis";
  if ( camp == 2) sim_type = "Gibbs";
  Input("input.dat",camp,temperature); //Inizialization

// EQUILIBRIZZAZIONE
  for(int i = 0; i < 1000; i++){
  Move(metro);
  ConfFinal();
}


  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
     Reset(iblk);   //Reset block averages
     for(int istep=1; istep <= nstep; ++istep)
     {
       Move(metro);   //metropolis or Gibbs algorithm
       Measure();    //Measure properties
       Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  tempGraf(sim_type, temperature);
  return 0;

}


void Input(string new_input,int camp,double temperature)
{
  ifstream ReadInput,Seed,ReadConf;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   //Read input informations
  ReadInput.open(new_input);

  ReadInput >> restart;		// se c'e 0 riparte dalla configurazione iniziale altrimenti no

  if(restart) Seed.open("seed.out");
  else Seed.open("seed.in");
  Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  Seed.close();
  
//Read input informations

  ReadInput >> temp;
  temp = temperature;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  metro = camp;       //non legge piu il campionamento dal file di input!!

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  cout << "Read initial configuration" << endl << endl;
  if(restart){
      ReadConf.open("config.final");
      for(int i = 0; i < nspin; i++)  ReadConf >> s[i];
  }
  else{
      for (int i=0; i<nspin; ++i)
      {
          if(rnd.Rannyu() >= 0.5) s[i] = 1;
          else s[i] = -1;
       }
  }
  ReadConf.close();
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)      //flippa( o non ) tutti gli spin presenti nel sistema
{
  int o,s_flip;
  double r;
  double p, energy_old, energy_new, sm , DeltaE,acc;
  double energy_up, energy_down;


  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);    //seleziona in modo random una particella

    if(metro==1) ///////////METROPOLIS
    {
      attempted += 1;
      //1. Modifico l'orientazione dello spin entratto
      s_flip = s[o] * (-1);
      //2. Calcolo DeltaE
      DeltaE = (-2.0)*Boltzmann((-1)*s_flip,o);
      //3.Probabilità che la mossa sia accettata
      p = exp((-1)*beta * DeltaE);
      acc = min<double>(1,p);
//cout<<acc<<endl;

      if ( acc == 1) {     //se l' energia è diminuita, accetta sempre lo scambio
        s[o] = s[o]*(-1);
        accepted +=1;
        } 
      else{
          //genera un numero random uniformemente distribuito
          r = rnd.Rannyu();;
          if(r <= p) {     //accetta la mossa di inversione dello spin
            s[o] = s[o]*(-1);
            accepted +=1;
          } 
      }
    }
    else ///////////////////GIBBS
    {
      //1.Estraggo un valore di spin random r ( +1 o -1)
      attempted += 1;
      r = rnd.Rannyu();
      if ( r <= 0.5) r = -1;
      if ( r > 0.5)  r = 1;
      //2. Calcolo DeltaE
      DeltaE = -2.0*Boltzmann(r,o);
      p = 1.0/(1.0+exp(-beta * DeltaE));
      if(rnd.Rannyu() < p){   //impongo il flip
        s[o] = r;
      }
      else  s[o] = -r;
      accepted += 1;
    }
  }
}

double Boltzmann(int sm, int ip)    //calcola il DeltaE*0.5 di energia a seguito di un flip. sm = iniziale dello spin che viene flippato; ip = posizione dello spin flippato
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0 ,u2 = 0.0, m = 0.0, c = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    u2 += u*u;
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = beta*(pow(m,2));
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
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
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
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

  // ENERGIA 
    Ene.open("output.ene.0",ios::app);    //blk_norm = numero di blocchi presi in considerazione fino a quel momento
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

  // CAPACITA TERMICA
    Heat.open("output.heat.0",ios::app);
    stima_c =beta*beta*( (blk_av[ic]/blk_norm) - pow(blk_av[iu]/blk_norm , 2)  )/(double)nspin;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

  //MAGNETIZZAZIONE
    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; 
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

  //SUSCETTIVITA MAGNETICA
    Chi.open("output.chi.0",ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; 
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

//cout << "Print final configuration to file config.final " << endl << endl;

  WriteConf.open("config.final");

  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void tempGraf(string sim_type,double temperature){  //calcola i valori medi relativi a tutti i blocchi, in funzione della temperatura
  ofstream Ene, Heat, Mag, Chi;
  const int wd=12;
//calcolo energia, calore specifico e suscettivita magnetica solo per h = 0!
  if ( h < 0.000001 ){
    // ENERGIA 
    Ene.open("output/" + sim_type+"/ene.dat",ios::app);
    Ene<< temperature<<setw(wd)<< glob_av[iu] /double(nblk)<<setw(wd)<<err_u<<endl;
    Ene.close();

  // CAPACITA TERMICA
    Heat.open("output/" + sim_type+"/heat.dat",ios::app);
    Heat << temperature <<  setw(wd) << glob_av[ic]/double(nblk)<< setw(wd) << err_c << endl;
    Heat.close();

  //SUSCETTIVITA MAGNETICA

    Chi.open("output/" + sim_type+"/chi.dat",ios::app);
    Chi << temperature <<  setw(wd) << glob_av[ix]/double(nblk)<< setw(wd) << err_x << endl;
    Chi.close();

  //MAGNETIZZAZIONE

    Mag.open("output/" + sim_type+"/mag_0.dat",ios::app);
    Mag << temperature <<  setw(wd) << glob_av[im]/double(nblk)<< setw(wd) << err_m << endl;
    Mag.close();

  }
  else{
//la magnetizzazione la calcolo in entrambe le casistiche!!
  //MAGNETIZZAZIONE
    Mag.open("output/" + sim_type+"/mag_2.dat",ios::app);
    Mag << temperature <<  setw(wd) << glob_av[im]/double(nblk)<< setw(wd) << err_m << endl;
    Mag.close();
  }

}


int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
