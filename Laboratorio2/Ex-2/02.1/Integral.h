#ifndef __INTEGRAL_H__
#define __INTEGRAL_H__

#include <algorithm>
#include <cstdlib>
#include <iomanip>

#include "random.h"
#include "FunzioniR1.h"

using namespace std;


		//Classe MONTECRALO contenente i metodi:

			//INTEGRAL_AVE: Metodo della media con distribuzione uniforme
			//INTEGRAL_AVE_IS: Calcolo dell' integrale mediante Importance Sampling


//////////////////////////////////////////////////////////classe che permette di integrare una funzione attraverso il METODO DI MONTECARLO
class Montecarlo{
public:
	Montecarlo(){};
	~Montecarlo(){};

	double IntegralAVE(double xmin,double xmax, unsigned int N, FunzioneBase*f,double partenza, Random *rand){	//calcola l' integrale con il metodo della MEDIA nell' intervallo [xmin, xmax] della funzione f, attraverso N punti
		
		double sum = partenza;
		double x;
		for(int i = 0; i <= N; i++){
			x = rand->Rannyu(xmin,xmax);
			sum += f->Eval(x);
		}
		return sum*(xmax-xmin)/(double)N;	//valore approssimato dell' integale
	} 

	double IntegralAVE_IS(double xmin,double xmax, unsigned int N,double partenza, Random *rand){	//calcola l' integrale con il metodo della MEDIA nell' intervallo [xmin, xmax] della funzione f, attraverso N punti
		
		double sum = partenza;
		double x;
		double num,den;
		for(int i = 0; i <= N; i++){
			x = rand->Linear(xmin,xmax);
			num = M_PI * 0.5 * cos(M_PI*0.5*x);
			den = 2 * (1-x);
			sum += (num/den);
			}
		return sum/(double)N;	//valore approssimato dell' integale
	} 

};

#endif


