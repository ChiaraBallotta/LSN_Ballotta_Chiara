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

//Simuliamo i numeri casuali
	int M = 10000;		//Total number of throws
	int N = 100;		//Total number of blocks
	double L = M/N;	//Number of throws in  each block

	double rand[M];		//Vettore che conterra i numeri casuali
	for ( int i=0; i<M; i++){ //creazione dei numeri casuali
	rand[i] = rnd.Rannyu();
	}
	
	//rnd.SaveSeed();		//salva in output il seme a cui siamo arrivati!
//calcoliamo i valori medi dei vari blocchi e li mandiamo su un file
	vector<double> ave;		//vettore che contiene la media di ogni blocco
	vector<double> ave2;

	int cont=0;		//segnatura della posizione nell'array
	for( int i=0;i<N;i++){	// i = indicatore del blocco n-esimo in cui ci troviamo
		double sum=0;	//contatore della somma per ciascun blocco
		for(int j=0;j<L;j++){
			cont= j+ i*L;
			sum = sum + rand[cont]; 
		}
		ave.push_back(sum/L);
		ave2.push_back( pow(sum/L,2));
	}
	//Print<double>("medie.dat",ave);

//calcolo la deviazione standard dei valori sul vector ave
	vector<double> dev;		// deviazione standard
	vector<double> dev2;
	vector<double> err; 		//incertezza statistica sulla media
cd 		double sum_dev2 = 0;
		for(int j = 0; j<i +1;j++){
			sum = sum + ave[j];
			sum2 = sum +  ave2[j];
			//cout<<"    j: "<<j<<"  "<<sum<<"  "<<sum2<<endl;
		}
		sum = pow(sum/i+1,2);
		//cout<<sum<<endl;
		sum2 = sum2/i+1;
		//cout<<sum2<<endl;
		dev.push_back(sqrt(abs(sum2-sum))) ;
		//cout<<"dev:  "<<dev[i-1]<<endl;
		sum_dev = sum_dev + dev[i];
		
		dev2.push_back(abs(sum2-sum));
		
		err.push_back(dev[i-1]/sqrt(i));
		err_dev.push_back()
		}
// stampo l'errore statistico e la media su un file
	vector<double> Num;
	for(int i = 1; i < N+1;i++){
		Num.push_back(L*i);
	}
	Print<double>("Medie.dat",Num,ave,err);


	//for(int i=0; i<20; i++){
      	//cout << rnd.Rannyu() << endl;
   //}

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
