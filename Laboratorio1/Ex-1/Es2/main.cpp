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

//////////////////// ESERCITAZIONE NUMERO 1.2


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

    int M = 10000;
	int N[4] = {1,2,10,100};

	//Inizializzo i file su cui salvare i risultati
	ofstream Unif,Exp,Lor;
	Unif.open("Unif.dat");
	Exp.open("Exp.dat");
	Lor.open("Lor.dat");
	//Ciclo sul numero di "Dadi"
	for(int i = 0; i < 4; i++){
		double summ_unif = 0;
		double summ_exp = 0;
		double summ_lor = 0;
		
		for(int j = 0; j < M; j++){		// ciclo sui lanci all' interno di ciscun dado
			summ_unif = 0;
			summ_exp = 0;
			summ_lor = 0;
			
			for(int k = 0; k < N[i]; k++){
				summ_unif += rnd.Rannyu();
				summ_exp += rnd.Exponential(1.);
				summ_lor += rnd.Cauchy(0.,1.);
			}
			summ_exp /= double(N[i]);
			summ_lor /= double(N[i]);
			summ_unif /= double(N[i]);

			Unif<<summ_unif<<endl;
			Lor<<summ_lor<<endl;
			Exp<<summ_exp<<endl;
		}
	}
	Unif.close();
	Exp.close();
	Lor.close();

    return 0;
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
