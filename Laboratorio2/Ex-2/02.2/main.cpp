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
#include <algorithm>
#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>

#include "random.h"
#include "RW.h"

using namespace std;


int main (int argc, char *argv[]){
//generatore numeri casuali
	Random rnd;  
    rnd.SetSeed();

	int M = 10000;		//Total number of throws
	int N = 100;		//Total number of blocks
	int L = M/N;		//Numero di passi in ciascun blocco

	int n = 100; 	//Numero di passi in ciascun RW
	int a = 1; 		//Passo del reticolo


//////////////////// ESERCIZIO 1 : Random Walks DISCRETI	
	Posizione_D *Random_Walk_D = new Posizione_D(rnd);	//Discreto
	Posizione_C*Random_Walk_C = new Posizione_C(rnd,a);	//Continuo
	
	// variabili utilizzate durante il Data-Blocking
	double Rm2_summ[n][N] = {0};	// Conterrà la distanza media, per passi fissati, relativa a ciascun blocco	
	double Rm2;
	int pos[]{0,0,0};
	double ave[n];
	double ave2[n];
	double dev = 0;

	double Rm2_summ_C[n][N] = {0};	// Conterrà la distanza media, per passi fissati, relativa a ciascun blocco	
	double Rm2_C;
	double pos_C[]{0,0,0};
	double ave_C[n];
	double ave2_C[n];
	double dev_C = 0;

// 1. Calcolo delle disastanze medie in ciascun blocco, per n fissato
	for( int i=0;i<N;i++){	// i = indicatore del blocco n-esimo in cui ci troviamo

		Posizione_D *Random_Walk_D = new Posizione_D(0,0,0,rnd);
		Posizione_C *Random_Walk_C = new Posizione_C(0,0,0,rnd,a);

		for(  int j = 0; j < L; j++){					//all' interno di ciascun blocco calcoliamo "L" RW, deve sempre ripartire da (0,0,0)
			
			for(int i = 0; i < 3 ; i++){
				pos[i] = 0;
				pos_C[i] = 0;
			} 

			for(int k = 0; k < n; k++){
				Random_Walk_D->Update(pos,rnd);
				Random_Walk_C->Update(pos_C,rnd);

				Posizione_D *New_Random_Walk_D = new Posizione_D(pos[0],pos[1],pos[2],rnd);
				Posizione_C *New_Random_Walk_C = new Posizione_C(pos_C[0],pos_C[1],pos_C[2],rnd,a);

				Rm2 = New_Random_Walk_D->Dist2();		//distanza del passo dall'origine
				Rm2_C = New_Random_Walk_C->Dist2();		//distanza del passo dall'origine
				// sommo la distanza a n fissato, per ciascun passo, per andare a mediarla sui blocchi
				Rm2_summ[k][i] += Rm2/double(L);
				Rm2_summ_C[k][i] += Rm2_C/double(L);
			}
		}
	}

	ofstream Discreto,Continuo;
	Discreto.open("RW_Discreto.dat");
	Continuo.open("RW_Continuo.dat");
// Calcolo della Medie nei blocchi

	for(int i = 0; i < n; i++){
		for(int j = 0 ; j < N; j++){
			double mean = Rm2_summ[i][j];
			double mean_C = Rm2_summ_C[i][j];

			ave[i] += mean/double(N);
			ave2[i] += mean*mean/double(N);
			ave_C[i] += mean_C/double(N);
			ave2_C[i] += mean_C*mean_C/double(N);
		}
	
		dev = sqrt((ave2[i])- pow(ave[i],2));
		dev_C = sqrt((ave2_C[i])- pow(ave_C[i],2));
		double err = 0;
		double err_C = 0;
		if(i == 0){	
			Discreto<<sqrt(ave[i])<<"  "<<0<<endl;
			Continuo<<sqrt(ave_C[i])<<"  "<<0<<endl;
		}
		else{
			err = dev/sqrt(N-1);
			Discreto<<sqrt(ave[i])<<"  "<<err<<endl;
			err_C = dev_C/sqrt(N-1);
			Continuo<<sqrt(ave_C[i])<<"  "<<err_C<<endl;
		}
	}
	Discreto.close();
	Continuo.close();
   //rnd.SaveSeed();		//salva in output il seme a cui siamo arrivati!!
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
