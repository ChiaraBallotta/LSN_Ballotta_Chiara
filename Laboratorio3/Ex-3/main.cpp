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

//////////////////// ESERCITAZIONE NUMERO 3 

template <typename T> void Print(const char * ,vector<T> ,vector<T>,vector<T>);


int main (int argc, char *argv[]){
//generatore numeri casuali
	Random rnd;
	rnd.SetSeed();	//inizializziamo il generatore dei numeri casuali

//Parametri
	double So = 100;		//prezzo asset a t = 0
	double T = 1;			//expire time
	double k = 100;			//strike price
	double r = 0.1;			//risk-free interest rate
	double sigma = 0.25;	// volatilità

//preparo la costruzione dei blocchi
	double M = 10000;		//campionamento totale
	double N = 100;		// numero dei blocchi
	double n = M/N;		//numeri di elementi per blocco

//costruzione contatore tempo, ESERCIZIO 3.2
	double time_steps = 100;
	double t = T/time_steps;		//distanza temporale tra un' intervallo e l'altro, sapendo che l'interavalo [0,T] è stato diviso in 100!!

	vector<double> Nf;

//////////////////ESERCIZIO 3.1:	CAMPIONAMENTO DIRETTO 
	double summC_block;
	double summP_block;
	double w;
	double s;
	double c;
	double p;
	double C_block;
	double P_block;

	double summC = 0;
	double summP = 0;
	double summC_2 = 0;
	double summP_2 = 0;
	double dev_C = 0;
	double dev_P = 0;

	vector<double> ave_C;		//valore medio in funzione dei blocchi del CALL
	vector<double> ave_P;		//valore medio in funzione dei blocchi del PUT
	vector<double> err_C;
	vector<double> err_P;	

//////////////////ESERCIZIO 3.2:	CAMPIONAMENTO DISCRETO
	double summC_block_d;
	double summP_block_d;
	double w_d;
	double s_d;
	double c_d;
	double p_d;
	double C_block_d;
	double P_block_d;

	double summC_d = 0;
	double summP_d = 0;
	double summC_2_d = 0;
	double summP_2_d = 0;
	double dev_C_d = 0;
	double dev_P_d = 0;

	vector<double> ave_C_d;		//valore medio in funzione dei blocchi del CALL
	vector<double> ave_P_d;		//valore medio in funzione dei blocchi del PUT
	vector<double> err_C_d;
	vector<double> err_P_d;

////////////  INIZIO SIMULAZIONE

	for( int j = 1; j < N + 1; j++){
		summC_block = 0;
		summP_block = 0;

		summC_block_d = 0;
		summP_block_d = 0;

		for ( int i = 1; i < (n + 1);i++){
		// Diretto
			w = rnd.Gauss(0,T);									
			s = So * exp(((r-sigma*sigma*0.5)*T) + sigma*w*pow(T,0.5));
			c = exp(-r*T) * max(0. ,s-k);		//call
			p = exp(-r*T) * max (0.,k-s);		//put
			summC_block = summC_block + c;
			summP_block = summP_block + p;

		// Discreto
			s_d = So;								//resetta il valore temporale inziale per il processo di discretizzazione
			for(int l = 0; l < time_steps; l++){
				w_d = rnd.Gauss(0,T);											
				s_d = s_d * exp(((r-sigma*sigma*0.5)*(t) + sigma * w_d * pow(t,0.5)));
			}
			c_d = exp(-r*T) * max(0. ,s_d-k);			//call
			p_d = exp(-r*T) * max (0.,k-s_d);			//put
			summC_block_d = summC_block_d + c_d;
			summP_block_d = summP_block_d + p_d;
		}

		C_block = summC_block/ double(n);		//calcola il valore medio all'interno di ciascun blocco CALL
		P_block = summP_block/double(n);		//calcola il valore medio all'interno di ciascun blocco PUT

		C_block_d = summC_block_d/ double(100);	// 100 = numero di intervalli in cui è stato diviso il periodo	
		P_block_d = summP_block_d/double(100);		

	//calcolo il valore medio in relazione al numero di blocchi	
		summC = summC + C_block;		//calcola la somma dei valori medi dei vari blocchi
		summP = summP + P_block;
		summC_2 = summC_2 + pow( C_block,2);		//calcola la somma dei valori medi dei vari blocchi
		summP_2 = summP_2 + pow(P_block,2);

		summC_d = summC_d + C_block_d;		
		summP_d = summP_d + P_block_d;
		summC_2_d = summC_2_d + pow( C_block_d,2);		
		summP_2_d = summP_2_d + pow(P_block_d,2);

		ave_C.push_back(summC /double(j));			//valore medio in funzione dei blocchi CALL
		ave_P.push_back(summP/double(j));			//valore medio in funzione dei blocchi PUT

		ave_C_d.push_back(summC_d /double(j));			
		ave_P_d.push_back(summP_d/double(j));		


		dev_C = sqrt((summC_2 /double (j)) - (pow(summC/double (j),2)));
		dev_P = sqrt((summP_2 /double (j)) - (pow(summP/double (j),2)));

		dev_C_d = sqrt((summC_2_d /double (j)) - (pow(summC_d/double (j),2)));
		dev_P_d = sqrt((summP_2_d /double (j)) - (pow(summP_d/double (j),2)));

		// Calcolo incertezza statistica
		if(j == 1) {
			err_C.push_back(0);
			err_P.push_back(0);
			}
		else{
		err_C.push_back(dev_C/sqrt(j-1));
		err_P.push_back(dev_P/sqrt(j-1));
		}

		if(j == 1) {
			err_C_d.push_back(0);
			err_P_d.push_back(0);
			}
		else{
		err_C_d.push_back(dev_C_d/sqrt(j-1));
		err_P_d.push_back(dev_P_d/sqrt(j-1));
		}

		Nf.push_back(j);
	}
	//Risultati ottenuti, stampati su file output
	Print<double>("Risultati_C.dat",Nf,ave_C,err_C);
	Print<double>("Risultati_P.dat",Nf,ave_P,err_P);

	Print<double>("Risultati_C_d.dat",Nf,ave_C_d,err_C_d);
	Print<double>("Risultati_P_d.dat",Nf,ave_P_d,err_P_d);
	//rnd.SaveSeed();		//salva in output il seme a cui siamo arrivati!

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


	





