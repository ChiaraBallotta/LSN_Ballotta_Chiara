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

#include "popolazione.h"
#include "random.h"
#include "input.h"


using namespace std;

template <typename T> void Print(const char * ,vector<T> ,vector<T>,vector<T>);
void cities_circ( int,Random,double**);	//genera la distribuzione delle citta
void cities_quadr(int, Random,double**);



int main (int argc, char *argv[]){

	int distr = 0;
	cout<<"Problema del commesso viaggiatore."<<endl;
	cout<<"Distribuzione uniforme delle città :"<<endl;
	cout<<"Selezionare 1 per distribuzione CIRCOLARE"<<endl<<"Selezionare 2 per distribuzione QUADRATICA:"<<endl;
	cin >> distr;
	while(distr != 1 && distr != 2){
		cout<<"Errore! Inserire un valore corretto"<<endl;
		cin >> distr;
	}
//generatore numeri casuali
	Random rnd;
	rnd.SetSeed();

//Generatore della distribuzione di città
	int rows = 3;
	int cols = N_cities;
	double** cities = new double*[rows];		// matrice che conterrà le città, e le cordinate x, y
	for (int i=0; i<rows; i++){
		cities[i] = new double[cols];
	}


	string type;

	if ( distr == 1){
		type ="Circolare";
		cities_circ(N_cities,rnd,cities);
	}
	if ( distr == 2){
		type ="Quadrato";
		cities_quadr(N_cities,rnd,cities);
	}
// Genero la popolazione
	Popolazione popolaz( N_individui, N_cities,rnd, cities );;
	popolaz.sorting();


// Inizio il processo di Mutazione, mediante N_generaz generazioni
	ofstream Best,Generaz,Best_ave,Best_best,Best_second;
	int ind = 0;

	Best.open(type + "/best.dat");			//stampa il migliore ogni generazione
	Generaz.open(type + "/generaz.dat");	//tiene traccia del progresso dell'algoritmo
	Best_ave.open(type + "/best_ave.dat");	// valore medio della lunghezza sulla prima migliore meta della popolazione
	Best_best.open(type + "/best_best.dat");	//contiene il migliore percrso in assoluto
	Best_second.open(type + "/best_second.dat");	//contiene il secondo migliore percrso in assoluto

	int summ;
	double ave;

	//Variabili per trovare il migliore e il secondo migliore
	double first[N_cities +1];
	first[0] = 100;
	double second[N_cities +1];
	second[0] = 100;

	for( int i = 0; i < N_generaz;i++){
		cout<<"GENERAZIONE N° "<<i+1<<endl;
		Generaz<<" Inizio generazione n° "<<i+1<<endl;
		popolaz.Mutation(rnd,cities, Generaz);

		Best<<popolaz.Get_individuo(0).Get_length()<<"  ";
		// calcolo del valore medio
		summ = 0;	
		for( int j = 0; j < N_individui/2.0;j++)	summ += popolaz.Get_individuo(j).Get_length();
		ave = summ / (N_individui/2.0);
		Best_ave<<i<<" "<<ave<<endl;

		//calcolo il migliore percorso in generale
		if(popolaz.Get_individuo(0).Get_length() <= first[0]) {
			first[0] = popolaz.Get_individuo(0).Get_length();
			for(int f = 1; f <N_cities +1;f++)	first[f] = popolaz.Get_individuo(0).Get_cities(f-1);
		}
		if(popolaz.Get_individuo(0).Get_length() > first[0] && popolaz.Get_individuo(0).Get_length() < second[0] ) {
			second[0] = popolaz.Get_individuo(0).Get_length();
			for(int f = 1; f <N_cities +1;f++)	second[f] = popolaz.Get_individuo(0).Get_cities(f-1);
		}

		for(int j = 0; j < N_cities; j++){		// stampa la migliore citta ci sciascuna mutazione su un file Best
			ind = popolaz.Get_individuo(0).Get_cities(j);	//indica la citta che stiamo prendendo in considerazione
			Best<<"  "<<ind<<" ";
		}
		Best<<endl;
		Generaz<<" Fine generazione n° "<<i+1<<endl;

	}
	Best.close();
	Generaz.close();
	Best_ave.close();
	for(int k = 0; k < N_cities + 1; k++){
		Best_best<<first[k]<<endl;
		Best_second<<second[k]<<endl;
	}
	Best_best.close();
	Best_second.close();



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


void cities_circ( int N,Random rnd,double** cities){		// genera le citta distribuite in modo casuale su una circonferenza di raggio 1
	//int rows = 3;
	int cols = N;
	double rand;

	ofstream out;
	out.open("Circolare/cities.dat");

	for(int i = 0; i < cols ; i++){		//N.B. trascuro la possibilita  che  me le generi nello stesso punto
		rand = rnd.Rannyu(0,2*M_PI);
//cout<<rand<<endl;
		cities[0][i] = i+1;
		cities[1][i] = cos(rand);		//cordinata x;
		cities[2][i] = sin(rand);		//cordinata y;

		out<< i+1<<"  "<<cities[1][i]<<"  "<<cities[2][i]<<endl;
	}
	out.close();	
}

void cities_quadr( int N,Random rnd,double** cities){		// genera le citta distribuite in modo casuale su una circonferenza di raggio 1
	//int rows = 3;
	int cols = N;
	double x,y;

	ofstream out;
	out.open("Quadrato/cities.dat");

	for(int i = 0; i < cols ; i++){		//N.B. trascuro la possibilita  che  me le generi nello stesso punto
		x = rnd.Rannyu();
		y = rnd.Rannyu();
		cities[0][i] = i+1;
		cities[1][i] = x;		//cordinata x;
		cities[2][i] = y;		//cordinata y;

		out<< i+1<<"  "<<cities[1][i]<<"  "<<cities[2][i]<<endl;
	}
	out.close();	
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
