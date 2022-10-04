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

//Programmazione in parallelo
#include "mpi.h"


using namespace std;
//using namespace arma;

template <typename T> void Print(const char * ,vector<T> ,vector<T>,vector<T>);
void cities_USA( double**);	//genera la distribuzione delle citta

void Salva_best_indv(ofstream& , ofstream&,ofstream&  , Popolazione ,double* ,int ,int,int );




int main (int argc, char *argv[]){

		
//generatore numeri casuali
	Random rnd;
	rnd.SetSeed();

// Scelta del tipo di Simmulazione
	cout<<"Problema del commesso viaggiatore."<<endl;
	cout<<" 50 Capitali degli Stati Uniti."<<endl;

	if(distr == 1){
		cout<<"Simulazione svolta mediante 1 continente"<<endl;
		cout<<"Assicurarsi di aver selezionato nel makefile che verrà utilizzato solo 1 nodo !!"<<endl;
	}
	if(distr == 2){
		cout<<"Simulazione svolta mediante Continenti Indipendenti"<<endl;
	}
	if(distr == 3){
		cout<<"Simulazione svolta mediante Continenti con Migrazione"<<endl;
	}

	//Generatore della distribuzione di città
	int rows = 3;
	int cols = N_cities;
	double** cities = new double*[rows];		// matrice che conterrà le città, e le cordinate x, y
	for (int i=0; i<rows; i++){
		cities[i] = new double[cols];
	}

	cities_USA(cities);

	// Genero la popolazione
	Popolazione popolaz( N_individui, N_cities,rnd, cities );
	popolaz.sorting();

	//Inizializzazione di MPI
	int size,rank;
	MPI_Init(&argc,&argv);					//Inizialisation of MPI enviroment
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status stat;

	//Inizializzo ciascun nodo su seed differenti
	rnd.Set_Primes(rank);

// Inizio il processo di Mutazione, mediante N_generaz generazioni
	ofstream Best,Best_ave,Best_best,Generaz,Nodi;

	Best.open(to_string(distr) + "/"+to_string(rank) +  "_best.dat");			//stampa il migliore ogni generazione
	Generaz.open(to_string(distr) + "/"+to_string(rank) +   "_generaz.dat");	//tiene traccia del progresso dell'algoritmo
	Best_ave.open( to_string(distr) + "/"+to_string(rank) +  "_best_ave.dat");	// valore medio della lunghezza sulla prima migliore meta della popolazione
	Best_best.open(to_string(distr) +"/"+to_string(rank) +   "_best_best.dat");	//contiene il migliore percrso in assoluto

	Nodi.open(to_string(rank) +"_nodi.dat");

	//Variabili per trovare il migliore e il secondo migliore
	double* first = new double[N_cities + 1];
	first[0] = 10000;

	int N_migraz = 10;    //numero di generazioni tra una mgrazione e l'altra
    int rep = N_generaz/N_migraz;           //numero di volte in cui il processo mutazioni + migrazione deve avvenire

	int contatore_mutaz = 0;					// conta il numero totali di mutazioni/migrazioni che avvengono

    for(int k = 0;  k < rep;k++){
        // 1. Mutazione
        for(int i = 0; i < N_migraz -1;i++){
            popolaz.Mutation(rnd,cities, Generaz);
			Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui,contatore_mutaz);
			contatore_mutaz += 1;
		}
		MPI_Barrier(MPI_COMM_WORLD);
Nodi<<" Esecuzione di "<<(k+1)*N_migraz<<" Migrazioni"<<endl;

        // 2. Migrazione
        if(distr != 3){
            popolaz.Mutation(rnd,cities,Generaz);
			Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui,contatore_mutaz);
			contatore_mutaz += 1;
        }
        else{
            int send = 0; 				// seleziona in modo Ranndom i due core in cui potrebbe avvenire la mutazione
			int receiver = 1;
			int scambio1[N_cities];		// conterrà l' individuo da scambiare	
			int itag = 1;
			int j = 0;  

            if( rank == 0){
				int p = 20;					//valore da me selezionato
				j = int(N_individui * pow(rnd.Rannyu(),p)) + 1;	//predilige i numeri piu bassi
				send = int(rnd.Rannyu(0,size));
				receiver = int(rnd.Rannyu(0,size));
				while( send == receiver)	receiver = int(rnd.Rannyu(0,size));
			}
			MPI_Barrier(MPI_COMM_WORLD);

            //Comunico a tutti nodi, i valori di j,send, receiver calcolati dal nodo zero
			MPI_Bcast(&j,1,MPI_INTEGER,0,MPI_COMM_WORLD);
			MPI_Bcast(&send,1,MPI_INTEGER,0,MPI_COMM_WORLD);
			MPI_Bcast(&receiver,1,MPI_INTEGER,0,MPI_COMM_WORLD);
Nodi<<" Migrazione numero : "<<k+1<< "  (j,send,receiver): "<<j<<" "<<send<<" "<<receiver<<endl; 

            if(rank != send && rank != receiver){
            	popolaz.Mutation(rnd,cities,Generaz);
				Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui,contatore_mutaz);
				contatore_mutaz +=1;
            }
			MPI_Barrier(MPI_COMM_WORLD);

            if(rank == send){
				for(int l = 0; l < N_cities; l++)	scambio1[l] = popolaz.Get_individuo(j).Get_cities(l);
				MPI_Send(scambio1,N_cities,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD);
            	popolaz.Mutation(rnd,cities,Generaz);
				Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui,contatore_mutaz);
				contatore_mutaz += 1;
			}			
			MPI_Barrier(MPI_COMM_WORLD);

			if(rank == receiver){
				MPI_Recv(scambio1,N_cities,MPI_INTEGER,send,itag,MPI_COMM_WORLD,&stat);
				//Elimino l' individuo peggiore della popolazione e salvo quello nuovo
				Individuo new_ind(scambio1,N_cities);
				new_ind.Set_length(cities);
				popolaz.Set_individuo(new_ind,popolaz.Get_dim() -1);
				popolaz.sorting();
				Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui,contatore_mutaz);
				contatore_mutaz += 1;
			}
			MPI_Barrier(MPI_COMM_WORLD);

			popolaz.sorting();	

			//Salva_best_indv(Best, Best_ave, Generaz, popolaz, first,N_cities,N_individui);
		}
Nodi<<"Il processo di migrazione si è concluso"<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	Best.close();
	Generaz.close();
	Best_ave.close();
	for(int k = 0; k < N_cities + 1; k++){
		Best_best<<first[k]<<endl;
	}
	Best_best.close();


	MPI_Finalize();

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


void cities_USA( double** cities){		// genera le citta distribuite in modo casuale su una circonferenza di raggio 1
	//int rows = 3;
	ifstream USA;
	USA.open("American_capitals.dat");
	if(!USA){
		cout<<"Cannot open file American_capitals.dat"<<endl;
		exit(11);
	}

	double val;
	int j = 0;

	while(!USA.eof()){
		cities[0][j] = j + 1;
		USA >> val;
		cities[1][j] = val;
		USA >> val;
		cities[2][j] = val;
		j += 1;
	}
	USA.close();	
}

// Funzione che salva il migliore di ogni generazione sui file di Interesse
void Salva_best_indv(ofstream& Best, ofstream& Best_ave, ofstream& Generaz, Popolazione popolaz,double* first,int N_cities,int N_individui,int contatore_mutaz){

		Best<<popolaz.Get_individuo(0).Get_length()<<"  ";
		// calcolo del valore medio
		int summ = 0;
		double ave;	
		int ind = 0;
		for( int j = 0; j < N_individui/2.0;j++)	summ += popolaz.Get_individuo(j).Get_length();
		ave = summ / (N_individui/2.0);
		Best_ave<<contatore_mutaz<<" "<<ave<<endl;

		//calcolo il migliore percorso in generale
		if(popolaz.Get_individuo(0).Get_length() <= first[0]) {
			first[0] = popolaz.Get_individuo(0).Get_length();
			for(int f = 1; f <N_cities +1;f++)	first[f] = popolaz.Get_individuo(0).Get_cities(f-1);
		}

		for(int j = 0; j < N_cities; j++){		// stampa la migliore citta ci sciascuna mutazione su un file Best
			ind = popolaz.Get_individuo(0).Get_cities(j);	//indica la citta che stiamo prendendo in considerazione
			Best<<"  "<<ind<<" ";
		}
		Best<<endl;
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
