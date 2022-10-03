#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Ex-8.h"

using namespace std;

int main(){  

    Input();

    x_old = x_i;

///////////////////////// ESERCIZIO 8.1  //////////////////////// 
    if(SA == 0) {
      Averages(1); //Calcolo di <H> su ogni blocco!
    }
    
///////////////////////// ESERCIZIO 8.2  //////////////////////// 
    if(SA != 0){

      ////// TEMPERATURA DI EQUILIBRIO--> MU E SIGMA MINIMIZZATI
      ofstream H_fin,Parametri,Text;
      Parametri.open("8.2/parametri.out");      // Valori di sigma e mu, in funzione dell'energia
      Text.open("8.2/text.out");                 //Descrizione di quello ch accade nei blocchi
      H_fin.open("8.2/energy.out");             // Energia in funzione della temperatura
      int blocco = 0;                         // indicatore del numero di  vvolte in cui viene abbassata la temperatura
      
      //EQUILIBRAZIONE PER AVVICINARMI AL PUNTO DI MINIMO
      for( int j = 0; j<100;j++) {
        Move_MC_metro();    // x_old --> punto di minimo dell'energia
      }
      cout<<x_old<<endl; 

//      while(err_H > deltaH ){      
        for( int l = 0; l < 100;l++){  
        
        double H_old, H_new, p,q;
        double diff_H = 0;

        H_old = globH_ave/double(nblk);

        for(int i = 0; i < stepSA; i++){  // per ogni temperatura, faccio metropoilis 

          // Calcolo del nuvo valore di H mediante un ciclo Montecarlo. N.B. bisogna cambiare i valori di mu  ìe sigma
          d_mu = rnd.Rannyu(-span_mu,span_mu)/beta;           // divido per beta siccome voglio che al diminuire della temperature, diminuisca anche di quanto varia
          d_sigma = rnd.Rannyu(-span_sigma,span_sigma)/beta;

          mu += d_mu;
          sigma += d_sigma;

          Averages(2);     // calcola l' energia in funzione dei blocchi, con i nuovi valori di sigma e mu

          H_new = globH_ave/double(nblk);

          // METROPOLIS      
          diff_H = H_new - H_old;
          p = exp(-beta*diff_H);
          q = min<double>(1,p);

          if ( diff_H <= 0){
            H_old = H_new;
          }
          else{
            if(rnd.Rannyu() < q){
               H_old = H_new;
            }
            else{   // No scambio: parametri ripristinati !!
              mu = mu - d_mu;
              sigma = sigma - d_sigma;
            }
          }
        }

        blocco = blocco + 1;

        H_fin <<beta<<"  "<<H_old <<"  "<<err_H<<endl;

        Parametri <<H_old<<"  "<<beta<<"  "<<mu<<"  "<<sigma<<endl;
        Text <<"-------------Blocco "<<blocco<<"-esimo---------------"<<endl;
        Text <<"beta : "<<beta<<endl<<"<H> : "<<H_old<<endl<<"Incertezza <H> : "<<err_H<<endl;
        Text <<"mu : "<<mu<<endl;
        Text <<"sigma : "<<sigma<<endl<<endl;
        beta += db;
      }

      Parametri.close();      // Valori di sigma e mu, in funzione dell'energia
      Text.close();                 //Descrizione di quello ch accade nei blocchi
      H_fin.close(); 

      ////// VALORE DI <H> CON ERRORE, ALLA CONDIZIONE DI MINIMIZZAZIONE
      //for( int j = 0; j++; j < 100) Move_MC_metro();
      Averages(1);
      
       
    }

  
  return 0;

}


void Input()
{
  ifstream ReadInput,Seed,ReadConf;

  cout<< "Single quantum particle in 1D           " << endl;
  cout<< "Monte Carlo simulation with Metropolis algorithm            " << endl << endl;
  cout<< "Simulated Annealing (SA) algorithm      " << endl << endl;
  cout<< "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout<< "The program uses h =1 and m =1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   //Read input informations
    ReadInput.open("input.dat");

    ReadInput >> SA;
    ReadInput >> restart;
    Seed.open("seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    
    Seed.close();
  
//Read input informations

  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> stepSA;

  cout << "The program perform Metropolis moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput >> db;
  ReadInput >> temp_o;
  beta = 1.0/temp_o;
  ReadInput >> x_i;
  ReadInput >> x0;
  ReadInput >> xf;

  ReadInput >> deltaH;

  ReadInput >> mu;
  ReadInput >> sigma;
  ReadInput >> span_mu;    
  ReadInput >> span_sigma; 

  ReadInput.close();

//Prepare arrays for measurements
    iu = 0; //Energy

//initial configuration
    if(restart){                      //legge i valori di mu e sigma dall'output dello step di montecarlo precedente
      ReadConf.open("config.final");
      ReadConf >> mu;
      ReadConf >> sigma;
  }
    else{
        ReadInput >> mu;
        ReadInput >> sigma;
    }

    ReadInput >> d_mu;
    ReadInput >> d_sigma;

    ReadInput.close();
    ReadConf.close();
    
}

double Gaussiana_combinaz(double x){
    return exp(-(x-mu)*(x-mu)/(2*sigma*sigma))+exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
}

double Eval_H(double x){            // restituisce il valore dell' energia risolvendo schrodinger               
    double exp1 = exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
    double exp2 = exp(-(x+mu)*(x+mu)/(2*sigma*sigma));
    double V = pow(x,4)-(5.0/2.0)*pow(x,2);         //Valore potenziale
    
    double H_k = 0.5*pow(1.0/sigma,2)*(exp1*(1-pow((x-mu)/sigma,2)) + exp2*(1-pow((x +mu)/sigma,2)));   //Calcolo analitico della derivata
    H_k = H_k / Gaussiana_combinaz(x);
    return H_k + V;
}

void Move_MC_metro()      //Mossa di Montecarlo, mediante metropolis
{
    double p_new,p_old,loss,q,r;

    span = xf-x0;
  // Seleziona in modo  Random uno spostamento
    step_dx = rnd.Rannyu(-span,span);         
  // Definizione della nuova posizione
    x_new = x_old + step_dx;
    while ( abs(x_new) > xf){        //assicura di rimanere nell' intervallo
      step_dx = rnd.Rannyu(-span,span);         
      x_new = x_old + step_dx;
    }
    attempted = attempted + 1;

  // Calcolo della probabilità
    p_new = pow(Gaussiana_combinaz(x_new),2);
    p_old = pow(Gaussiana_combinaz(x_old),2);
    loss = p_new/p_old;

    q = min<double> ( 1,loss);

    if ( loss >= 1 ){  //il passo viene accettato
      x_old = x_new;
      accepted = accepted + 1;
    }

    else{
      r = rnd.Rannyu();
      if(r < q){
        x_old = x_new;
        accepted = accepted + 1;
      }
    }
}


void Averages(int es){ //Print results for current block. 
////////es = 1--> stampa in 8.1/8.2

/////// es = 2 --> non stampa energia e traiettoria

  ofstream H_ave, parameters,pos;
  int wd = 12;
  double H;
  string Ex;
  if(SA == 0)  Ex = "8.1";
  if ( SA == 1) Ex = "8.2";

// 1. <H> IN FUNZIONE DEL NUMERO DEI BLOCCHI, NELLA CONDIZIONE DI MINIMIZZAZIONE
  if(es == 1){
  H_ave.open(Ex +"/results.out");
  pos.open(Ex +"/pos.out");
  }

  globH_ave = 0;
  globH_ave2 = 0;

  for( int i = 0; i < nblk; i++){   // ciclo sul numero dei blocchi
    blkH_ave = 0;       // annullo il valore di <H>, all' inizio di ogni blocco
    attempted = 0;
    accepted = 0;

    if(es == 1) pos<<x_old<<setw(wd)<<Eval_H(x_old)<<endl;
    //inizio il ciclo all' interno di un singolo blocco per il calcolo di <H>
    for(int j = 0; j < nstep; j++){
      // Per ogni step fa una mossa Metropolis
      Move_MC_metro();
      // accumulo il valore dell' H calcolato per ogni step ( Measure() + Accumulate() in Ex-6)
      H = Eval_H(x_old);
      blkH_ave += H; 
      if(es == 1 && x_old == x_new) pos<<x_old<<endl;
    }
    
    // Calcolo medie sul blocco
    stima_H = blkH_ave / double(nstep);
    globH_ave +=stima_H;
    globH_ave2 +=pow(stima_H,2);
    if(i == 0)  {err_H = 0;}
    else{err_H = Error(globH_ave/(i+1),globH_ave2/(i+1),i+1);}

    if(es == 1) H_ave<< setw(wd)<<(i+1)<< setw(wd)<< globH_ave/double(i+1)<<setw(wd)<<err_H<<endl;
    
  }
  if(es == 1){
  pos.close();
  H_ave.close();  
  }
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}










