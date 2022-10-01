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

	int M = 10000;		//Total number of throws
	int N = 100;		//Total number of blocks
	int L = M/N;		//Number of throws in  each block
	vector<double> Nf;
	
///////// 1.1 Calcolo del Valore medio in ciascun Blocco
	double summ_block;
	double summ = 0;
	double summ2 = 0;
	double dev = 0;
	double ave;

	vector<double> ave_v;		//media per passo i-esimo, presi in considerazione i blocchi precedenti
	vector<double> err;

///////// 1.2 Valore medio della deviazione standard
	vector<double> ave_dev_v;
	vector<double> err_dev;
	double summ_dev_block;
	double summ_dev = 0;
	double summ2_dev = 0;
	double dev_dev = 0;
	double ave_dev;

	double r;
	for( int i = 1; i < N +1; i++){
		summ_block = 0;
		summ_dev_block = 0;
;
		for(int j = 0; j < M; j++){
			r = rnd.Rannyu();
		} 

		for(int j = 1; j < L +1 ; j++){
			r = rnd.Rannyu();
			summ_block = summ_block + r;
			summ_dev_block = summ_dev_block + pow(rnd.Rannyu()-0.5,2);
		}
		ave = summ_block/double(L);
		ave_dev = summ_dev_block/double(L);

		summ = summ + ave;
		summ2 = summ2 + pow(ave,2);
		summ_dev = summ_dev + ave_dev;
		summ2_dev = summ2_dev + pow(ave_dev,2);

		ave_v.push_back(summ/double(i));
		dev = sqrt((summ2/double(i))-pow(summ/double(i),2));
		ave_dev_v.push_back(summ_dev/double(i));
		dev_dev = sqrt((summ2_dev/double(i))-pow(summ_dev/double(i),2));
		//calcolo Incertezza statistica
		if(i == 1){
			err.push_back(0);
			err_dev.push_back(0);
		}
		else{
			err.push_back(dev/sqrt(i-1));
			err_dev.push_back(dev/sqrt(i-1));
		}
		Nf.push_back(i);
	}

///////// 1.3 Test del chi^2

	vector<double> chi2;
	vector<double> null;
	double chi = 0;
	M = pow(10,4);
	N = 100;
	double intervallo = 1.0/double(N);

	for(int i = 0; i < N;i++){
		chi = 0;
	    vector<double> occupazione;
		//Ciclo sui lanci in ogni intervallo
		for(int j = 0; j < M;j++){
			r = rnd.Rannyu();
			for( int k = 0; k< N;k++){              //Ciclo su tutti gli N blocchi
				occupazione.push_back(0);
				if(r <= 1.0/double(N*(k+1)) && r > 1.0/double(M*1) ){
					 occupazione[k] = occupazione[k] + 1;
				}
			}
		}
		for( int j = 0; j < N; j++){
		    chi += (pow((occupazione[j] - double(M)/double(N)),2))/(double(M)/double(N));
		}
		chi2.push_back(chi);
	    null.push_back(0);
	}

	//Stampo i risultati su file output
	Print<double>("Ave.dat",Nf,ave_v,err);
	Print<double>("Dev.dat",Nf,ave_dev_v,err_dev);
	Print<double>("Chi.dat",Nf,chi2,null);
	


    return 0;
}

template <typename T> void Print(const char * Filename,vector<T> data,vector<T> data2,vector<T> data3 ){
	 ofstream out;
	 out.open(Filename);
	 if(!out){
			cout<<"Non posso creare il file "<<Filename<<endl;
			return;
		}
		cout<<data.size()<<endl<<data2.size()<<endl<<data3.size()<<endl;
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
