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

#include "Integral.h"

using namespace std;

template <typename T> void Print(const char * ,vector<T> ,vector<T>,vector<T>);


int main (int argc, char *argv[]){
//generatore numeri casuali
	Random *rnd = new Random();  
    rnd->SetSeed();

//Simuliamo i numeri casuali
	int M = 10000;		//Total number of throws
	int N = 100;		//Total number of blocks
	double L = M/N;	//Number of throws in  each block

	double pi2 =  M_PI/double(2);
	FunzioneBase *f = new Coseno(pi2,pi2,0);
	Montecarlo *Integral = new Montecarlo();	//creazione della classe montecarlo
	double dx = 1;
	double sx = 0;
	int n_tiri = M;							//numero di tiri fatti all' interno di ciascun montecarlo

//////////////////// ESERCIZIO 1 : Calcolo il valore dell' integrale mediante distribuzione UNIFORME 

	vector<double> ave;		//vettore che contiene la media di ogni blocco con i blocchi precedenti!!
	vector<double> ave_dev;

	double Integral_MC;
	double summ_blocchi = 0;
	double summ_blocchi2 = 0;

	double dev2 = 0;
	double dev = 0;
	double err = 0;

//////////////////// ESERCIZIO 2 : Calcolo il valore dell' integrale mediante IMPORTANCE SAMPLING

	vector<double> ave_IS;		
	vector<double> ave_dev_IS;

	double Integral_MC_IS;
	double summ_blocchi_IS = 0;
	double summ_blocchi2_IS = 0;

	double dev2_IS = 0;
	double dev_IS = 0;
	double err_IS = 0;

	for( int i=0;i<N;i++){	// i = indicatore del blocco n-esimo in cui ci troviamo

		Integral_MC = Integral->IntegralAVE(sx,dx,n_tiri,f, 0,rnd);
		summ_blocchi += Integral_MC;
		summ_blocchi2 += pow(Integral_MC,2);

		Integral_MC_IS = Integral->IntegralAVE_IS(sx,dx,n_tiri, 0,rnd);
		cout<<Integral_MC_IS<<endl;
		summ_blocchi_IS += Integral_MC_IS;
		summ_blocchi2_IS += pow(Integral_MC_IS,2);

		ave.push_back(summ_blocchi / double(i+1));
		ave_IS.push_back(summ_blocchi_IS / double(i+1));

		dev = summ_blocchi / double(i+1);
		dev2 = summ_blocchi2/double(i+1);
		if(i == 0) {err = 0;}
		else {err = sqrt((dev2 - dev*dev) /double(i));}
		ave_dev.push_back(err);
	
		dev_IS = summ_blocchi_IS / double(i+1);
		dev2_IS = summ_blocchi2_IS/double(i+1);
		if(i == 0) {err_IS = 0;}
		else {err_IS = sqrt((dev2_IS - dev_IS*dev_IS) /double(i));}
		ave_dev_IS.push_back(err_IS);
	}
	
// stampo l'errore statistico e la media su un file
	vector<double> Num;
	for(int i = 1; i < N+1;i++){
		Num.push_back(i);
	}
	Print<double>("Medie.dat",Num,ave,ave_dev);
	Print<double>("Medie_IS.dat",Num,ave_IS,ave_dev_IS);


   //rnd.SaveSeed();		//salva in output il seme a cui siamo arrivati!!
   return 0;
}


template <typename T> void Print(const char * Filename,vector<T> data,vector<T> data2,vector<T> data3 ){
	 ofstream out;
	 out.open(Filename);
	 if(!out){
			cout<<"Non posso creare il file "<<Filename<<endl;
			return;
		}
		//cout<<data.size()<<endl<<data2.size()<<endl<<data3.size()<<endl;
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
