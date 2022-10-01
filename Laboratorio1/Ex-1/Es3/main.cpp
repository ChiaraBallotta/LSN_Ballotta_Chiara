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
#include <string>
#include <cstdlib>
#include <vector>
#include "random.h"
#include <algorithm>
#include <numeric>
#include <cmath>

using namespace std;

//////////////////// ESERCITAZIONE NUMERO 1

template <typename T> void Print(const char * ,vector<T> ,vector<T>,vector<T>);


int main (int argc, char *argv[]){

//generatore numeri casuali
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
        Primes >> p1 >> p2 ;
        } else cerr << "PROBLEM: Unable to open Primes" << endl;
        Primes.close();

        ifstream input("seed.in");
        string property;
        if (input.is_open()){
        while ( !input.eof() ){
        	input >> property;
        	if( property == "RANDOMSEED" ){
            		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            		rnd.SetRandom(seed,p1,p2);
         	}
     	}
        input.close();
        } else cerr << "PROBLEM: Unable to open seed.in" << endl;

//Simulazione dei Numeri Casuali

	int N = 100;		//Total number of blocks
	int M = 50000;		//Number of throws in  each block
	vector<double> Nf;
	
// costruzione del sistema
	double L = 3.0;
	double d = 5.0;
	int N_hit = 0;		// contatore del numero di volte in cui la barretta tocca gli estremi
	double y; 			//posizione centrale barretta
	double sin_theta;		// orientamento barretta
	double y_dx,y_sx;	//posizione estremi barretta

// Numeri casuali usati
	double a,b;

///////// Calcolo del Valore medio in ciascun Blocco
	double summ = 0;
	double summ2 = 0;
	double dev = 0;
	double ave;

	vector<double> ave_v;		//media per passo i-esimo, presi in considerazione i blocchi precedenti
	vector<double> err;

	for( int i = 1; i < N +1; i++){
		N_hit = 0;

		// Ciclo sugli M lanci
		for(int j = 1; j < M +1 ; j++){
			y = rnd.Rannyu(0,d);

			// Generazione dell' angolo casuale
			a = rnd.Rannyu(0,L);
			b = rnd.Rannyu(0,L);
			while(b > sqrt(pow(L,2)-pow(a,2))){
				b = rnd.Rannyu(0,L);
			}
			// Costruzione posizione ago
			sin_theta = b/sqrt(pow(a,2)+pow(b,2));
			y_dx = y + (L)*sin_theta;
			//y_sx = y - (L/2.0)*sin_theta;
			//if(y_dx < y_sx){
				//double t = y_sx;
				//y_sx = y_dx;
				//y_dx = t;
			//}
			// Controlla se interseca la linea in y = 0
			if(y_dx >= d || y_dx <= 0){		// Condizione in cui interseca la linee
				N_hit += 1;
			}
		}
		ave = (2.0*L*double(M))/double(N_hit*d);
		summ = summ+ave;
		summ2 = summ2 + pow(ave,2);

		ave_v.push_back(summ/double(i));
		dev = sqrt((summ2/double(i))-pow(summ/double(i),2));
		//calcolo Incertezza statistica
		if(i == 1){
			err.push_back(0);
		}
		else{
			err.push_back(dev/sqrt(i-1));
		}
		Nf.push_back(i);
	}


	Print<double>("Pi.dat",Nf,ave_v,err);
    return 0;
}

template <typename T> void Print(const char * Filename,vector<T> data,vector<T> data2,vector<T> data3 ){
	 ofstream out;
	 out.open(Filename);
	 if(!out){
			cout<<"Non posso creare il file "<<Filename<<endl;
			return;
		}
	if (data.size() == data2.size() && data2.size() == data3.size()){
   		 for(int i = 0; i < data.size(); i ++) {
				out<<data[i]<<" "<<data2[i]<<" "<<data3[i]<<endl;
		}
	}
	else{
		cout<<"Errore, i 3 array hanno dimensione diversa!"<<endl;
		return;
	}
    out.close();
 
  return;
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
