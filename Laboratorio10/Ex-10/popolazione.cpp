#include <cmath>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cstdlib>
#define NDEBUG 
#include <assert.h>

#include "popolazione.h"

using namespace std;


vector<int> blocco1;
vector<int> blocco2;
vector<int> indiv;

int rand1;
int rand2;

//				a. DISTANZA TRA DUE PUNTI
double dist(double x1, double y1,double x2,double y2){	
	double d = 0;
	d = pow(x1 - x2, 2) + pow(y1 - y2,2);
	return sqrt(d);
}


//				b.CHECK SU INDIVIDUO ( una volta che è stato creato)---> FALSE = ERRATO!!
bool check(Individuo ind){	
	indiv.clear();										
	int summ = 0;
	for(int i = 0; i < ind.Get_n();i++){	// copio il vettore di ingresso in un nuovo vettore
		indiv.push_back(ind.Get_cities(i));
		summ += indiv[i];
	}
	vector<int> null;
	for(int i = 0; i < indiv.size();i++) null.push_back(0);

	double index = 0;
	if(indiv[0] =! 1) {
		cout<<" Il primo elemento è diverso da zero !!"<<endl;
		return false;
		}
	if(ind.Get_n() != N_cities){
		cout<<"Numero di citta sbagliate!!"<<endl;
		return false;
	}
	if( summ != (N_cities*(N_cities +1))/2.0){
		cout<<"La somma dei numeri non è corretta --> è presente una ripetizione !!"<<endl;
		return false;
	}

	return true;
}



//				c. SCAMBIA LA POSIZIONE DI DUE INDIVIDUI
void scambiaByRef(Individuo& a, Individuo& b) {
	Individuo t = a;
	a = b;
	b = t;
}



/////////////////////		INDIVIDUO		//////////////////

Individuo::Individuo(int* new_ind,int dim){
	for(int i = 0; i < dim;i++){
		m_ind.push_back(new_ind[i]);
	}
	m_n = dim;						//non soo sicura che sia corretta
}

//Genera un nuovo individuo in maniera Random
Individuo::Individuo(int n,Random& rnd){
	m_n = n;

	//riempiamo casualmente il vettore m_ind
	m_ind.push_back(1);	// il primo elemento deve sempre essere 1

	vector<int> sequenza;
	for(int i = 0; i < n-1; i++){
		sequenza.push_back(i+2);
	}
	
	int rd,cont;
	cont = n -2;
	for(int c = 1; c < n;c++){

		rd = rnd.Rannyu(0,cont);
		m_ind.push_back(sequenza[rd]);
		sequenza.erase(sequenza.begin() +rd);
		cont = cont-1;
	}
	//controllo che la dimensione di individuo sia corretta
	if(m_ind.size() != n) cout<<"Errore nella costruzione di individuo!";
}


//distruttore
Individuo::~Individuo(){
	m_ind.clear();
	m_ind.resize(0);
}

//Modifica un individuo, creandone uno nuovo Random
void Individuo :: Individuo_new(Random& rnd){

	m_ind[0] = 1;	// il primo elemento deve sempre essere 1

	vector<int> sequenza;
	for(int i = 0; i < this->Get_n() -1; i++){
		sequenza.push_back(i+2);
	}
	
	int rd,cont;
	cont = this->Get_n() -2;
	for(int c = 1; c < this->Get_n();c++){

		rd = rnd.Rannyu(0,cont +1);
		while ( rd ==  cont + 1)			rd = rnd.Rannyu(0,cont +1); 	// aggiunto per aumentare la probabilità di pescare l'ultimo
		m_ind[c] = sequenza[rd];
		sequenza.erase(sequenza.begin() +rd);
		cont = cont-1;
	}
	//controllo che la dimensione di individuo sia corretta
	if(m_ind.size() != this->Get_n()) cout<<"Errore nella costruzione di individuo!";
	return;
}



//metodi
void Individuo::Set_length(double** cities){
	double length = 0;
	int indx;
	double city1_x,city1_y,city2_x,city2_y;
	indx = m_ind[0];
	city2_x = cities[1][indx-1];
	city2_y = cities[2][indx-1];
		for( int i = 1; i < m_n ; i++){
			indx = m_ind[i];

			city1_x = city2_x;
			city1_y = city2_y;
			city2_x = cities[1][indx-1];
			city2_y = cities[2][indx-1];
			length = length + dist(city1_x,city1_y,city2_x,city2_y);
		}
	indx = m_ind[0];
	length = length + dist (city2_x,city2_y,cities[1][indx- 1],cities[2][indx -1]);
	m_length = length;
	return;
};

double Individuo :: Get_length(){ return m_length;}

int Individuo ::  Get_cities(int i){
	return m_ind[i];
}

int Individuo ::  Get_n() {return m_n;}

void Individuo :: Modif_city(int i ,int nuova){	//sostituisco la città i-esima con nuova
	m_ind[i] = nuova;
	return;
}

int Individuo :: Get_size(){ return m_ind.size();}
	





//////////////////		POPOLAZIONE		////////////////////
//generatore: crea una nuova popolazione, con individui definiti in maniera Random
Popolazione::Popolazione(int dim, int n,Random& rnd, double** cities ){

	m_dim = dim;
	ofstream Pop;
	Pop.open("popolaz_iniz.dat");		// file che contiene la popolazione iniziale

	Individuo indv(n,rnd);			//creo ogni volta un nuovo individuo in modo Random


	for(int i = 0; i< dim;i++){

		indv.Individuo_new(rnd);			//creo ogni volta un nuovo individuo in modo Random

		indv.Set_length(cities);

		m_popolaz.push_back(indv);
		m_popolaz_new.push_back(indv);

		for(int j = 0; j < n;j++)	Pop<<m_popolaz[i].Get_cities(j)<<"  ";
		Pop<<"lenght : "<<m_popolaz[i].Get_length()<<endl;

	}

	Pop.close();
}

Popolazione :: ~Popolazione(){};

//Inserisce l' individuo ind nella posizione i-esima
void Popolazione::Set_individuo(Individuo ind,int i){
		m_popolaz[i] = ind;
}

//ordina in ordine decrescente gli individuoi, in base alla length. Quelli con length piu picola sono alla alla fine
void Popolazione::sorting(){

	for(int i = 0; i < m_popolaz.size(); i++)
	for(int j = i + 1; j < m_popolaz.size(); j ++)
		if(m_popolaz[j].Get_length()<m_popolaz[i].Get_length())		
			scambiaByRef(m_popolaz[i],m_popolaz[j]);
}

int Popolazione :: Get_dim(){return m_dim;}


//Seleziona un individuo all'interno della popolazione. Da adoperare dopo aver fatto il sorting della popolazione
int Popolazione::selez(Random& rnd){		//restituisce la posizione dell'individuo selezionato
	int j;
	int p = 6;					//valore da me selezionato
	double rand = rnd.Rannyu();
	j = int(m_dim * pow(rand,p)) + 1;	//predilige i numeri piu bassi
	return j;
}

// Restituisce l' individuo i-esimo dal bach degli individui
  Individuo Popolazione ::  Get_individuo(int i){
		return m_popolaz[i];
  }


// 1) PERMUTA  due citta all'interno di un individuo, denominato dall'indice i
void Popolazione :: Pair_permutaz(int i,Random& rnd,double** cities){

	rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n() -1);					//indice della prima citta da scambiare 1
	rand2 = rnd.Rannyu(1,m_popolaz[i].Get_n()-1);					// indice della secondfa citta da scambiare 2
	while( rand2 == rand1) rand2 = rnd.Rannyu(1,m_popolaz[i].Get_n()-1);
	if ( rand1 >= m_popolaz[i].Get_n() ||rand2 >= m_popolaz[i].Get_n() ) cout<<"Problema!!"<<endl;
	double  scamb = m_popolaz[i].Get_cities(rand1);
	m_popolaz[i].Modif_city(rand1,m_popolaz[i].Get_cities(rand2));
	m_popolaz[i].Modif_city(rand2,scamb);
	
	m_popolaz[i].Set_length(cities);

	if(!check(m_popolaz[i])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};

	return;
}

// 3) PERMUTA UN BLOCCO DI n CITTA, all'interno di individuo, denominato dall'indice i-esimo
void Popolazione:: permutaz_N(int i, Random& rnd,double** cities){

	int n = rnd.Rannyu(1,(m_popolaz[i].Get_n()-1)/2.0);
	if(n >= (m_popolaz[i].Get_n()-1)/2.0)	cout<<"Errore! la permutazione_N non è eseguibile"<<endl;
	rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n()-1);			//indice inizio primo blocco
	while ((rand1 + n -1) >= (m_popolaz[i].Get_n() -1)) rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n() -1);
	//creo un vettore con tutti zero e 1 nelle caselle che subiranno la permutazione
	rand2 = rnd.Rannyu(1,m_popolaz[i].Get_n()-1);			//indice fine secondo blocco	
	while ( (rand1 + n -1) >= (m_popolaz[i].Get_n()-1) || abs(rand1 - rand2) < n) rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n());

	// Procedo coon la permutazione. Identifico i due blocchi che devo scambiare
	for(int j = 0; j < n;j++)	blocco1.push_back(m_popolaz[i].Get_cities(rand1+j));
	for(int j = 0; j < n;j++)	blocco2.push_back(m_popolaz[i].Get_cities(rand2+j));

	for(int j = 0; j < n; j++){
		m_popolaz[i].Modif_city(rand1 + j,blocco1[j]);
		m_popolaz[i].Modif_city(rand2 + j,blocco2[j]);

	}
	blocco1.clear();
	blocco2.clear();

	m_popolaz[i].Set_length(cities);

	if(!check(m_popolaz[i])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};	

	return;
}

// 2) SHIFTA m CITTA, DI n POSIZIONI; dell' inndividuo i-esimo
void Popolazione :: Shift(int i,Random& rnd,double** cities){

	int n = rnd.Rannyu(1,m_popolaz[i].Get_n()-2);
	int m;
	if( n > m_popolaz[i].Get_n()/2.0){
		m = rnd.Rannyu(1,m_popolaz[i].Get_n()-n-1);
	}
	else{
		m = rnd.Rannyu(1,n);
	}
	if (m >= (m_popolaz[i].Get_n()-1) - 1 )	cout<<"Errore! la permutazione Shift non è eseguibile"<<endl;
	rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n()-n-m);	

	// Procedo con lo shifting
	for(int j = 0; j < m;j++){
		blocco1.push_back(m_popolaz[i].Get_cities(rand1+j));
	}
	for(int j = 0; j < m;j++){
		blocco2.push_back(m_popolaz[i].Get_cities(rand1 +j + n));
		m_popolaz[i].Modif_city(rand1 +j + n,blocco1[j]);
		m_popolaz[i].Modif_city(rand1 + j,blocco2[j]);
	}
	blocco1.clear();
	blocco2.clear();

	m_popolaz[i].Set_length(cities);

	if(!check(m_popolaz[i])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};

	return;
}

// 4) INVERTE L'ORDINE DI M CITTA, dell' individuo i-esimo
void Popolazione:: Inversion(int i,Random& rnd,double** cities){

	rand1 = rnd.Rannyu(1,m_popolaz[i].Get_n());		// non puo essere l'ultima citta
	rand2 = rnd.Rannyu(rand1 + 1,m_popolaz[i].Get_n()-1);

	int scamb = 0;
	for( int j = rand1; j <= rand2;j++)	{
		blocco1.push_back(m_popolaz[i].Get_cities(j));
		scamb +=1;
		}
	for( int j = 0; j < scamb; j++){
		m_popolaz[i].Modif_city(rand1 + j,blocco1[scamb-1-j]);
	}
	blocco1.clear();

	m_popolaz[i].Set_length(cities);

	if(!check(m_popolaz[i])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};

	return;
}

// CROSSOVER tra l' individuo i-esimo e quello j-esimo
void Popolazione :: Crossover(int i, int j, Random& rnd,double** cities){

	rand1 = rnd.Rannyu(2,m_popolaz[i].Get_n()-2);		// indica il primo elemento che subira crossover. N.B. il cross over non puo partire dal secondo elemento altrimenti i figli sono uguali ai genitori, anaolgamento con l'ultimo
	vector<int> copy;

	for(int l = 0; l < m_popolaz[i].Get_n();l++){
		copy.push_back(m_popolaz[i].Get_cities(l));
	}
	for(int l = rand1; l < m_popolaz[i].Get_n(); l++){		// salvo i blocchi su fui fare crossover su vettori esterni
		blocco1.push_back(m_popolaz[i].Get_cities(l));
		blocco2.push_back(m_popolaz[j].Get_cities(l));
	}
	//adopero il crossover come indicato sul jupyter
	int selez1;
	int contatore1 = 0;
	for(int l = 1; l < m_popolaz[i].Get_n(); l++ ){				// modifico l' individuo i-esimo
		selez1 = m_popolaz[j].Get_cities(l);
		for(int k = 0; k < blocco1.size();k++){
			if( selez1 == blocco1[k]){
				m_popolaz[i].Modif_city(rand1 + contatore1,selez1);
				contatore1 += 1;
			}
		}
	}
	contatore1 = 0;
	for(int l = 1; l < m_popolaz[j].Get_n(); l++ ){				// modifico l' individuo j-esimo
		selez1 = copy[l];
		for(int k = 0; k < blocco2.size();k++){
			if( selez1 == blocco2[k]){
				m_popolaz[j].Modif_city(rand1 + contatore1,selez1);
				contatore1 += 1;
			}
		}
	}
	copy.clear();
	blocco1.clear();
	blocco2.clear();

	m_popolaz[i].Set_length(cities);
	m_popolaz[j].Set_length(cities);
	if(!check(m_popolaz[i])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};

	if(!check(m_popolaz[j])){
		cout<<"ERRORE: Individuo non costruito correttamente. Sistemare il Programma"<<endl;
	};

	return;
}


//Algoritmo che esegue le MUTAZIONI
void Popolazione :: Mutation(Random& rnd,double** cities,ofstream& Generaz){
	//Ordino la popolazione
	this->sorting();

	int ind1,ind2;
	double p,q1,q2;
	int contatore = 0;
	 while(contatore <= m_dim ){
		// 1) Selezioniamo due individui
		ind1 = this->selez(rnd);
Generaz<<"	mutazione "<<contatore+1<<endl;
		while(ind1 == 1 || ind1 == m_dim)	ind1 = this->selez(rnd);	// non voglio modificare il migliore, ovvero il primo

		ind2 = this->selez(rnd);
		while(ind2 == 1 || ind2 == m_dim)	ind2 = this->selez(rnd);

		// 1.1) Salvo i due individui su m_popolaz
		m_popolaz_new[contatore] = m_popolaz[ind1];
		contatore += 1;
		if(contatore == m_dim)	break;
		m_popolaz_new[contatore] = m_popolaz[ind2];
		contatore += 1;
		if(contatore == m_dim)	break;
		// 2) Opero Crossover sui due individui con una certa probabilità
Generaz<<"			Crossover"<<endl;
		p = rnd.Rannyu();

		if(p < probability[4] && ind1!=ind2){
				Crossover(ind1,ind2,rnd,cities);
				m_popolaz_new[contatore] = m_popolaz[ind1];
				contatore += 1;
				if(contatore == m_dim)	break;
				m_popolaz_new[contatore] = m_popolaz[ind2];
				contatore += 1;
				if(contatore == m_dim)	break;
		}

		// 3) Opero Mutazione sui figli ( se ne permuto uno, non è detto che permuto anche l'altro!!)
Generaz<<"			PairPermutation"<<endl;
		q1 = rnd.Rannyu();
		q2 = rnd.Rannyu();
		if(q1 < probability[0])	{
			this->Pair_permutaz(ind1,rnd,cities);
			m_popolaz_new[contatore] = m_popolaz[ind1];
			contatore += 1;
			if(contatore == m_dim)	break;
		}
		if(q2 < probability[0])	{
			this->Pair_permutaz(ind2,rnd,cities);
			m_popolaz_new[contatore] = m_popolaz[ind2];
			contatore += 1;
			if(contatore == m_dim)	break;
		}


Generaz<<"			Shift"<<endl;
		q1 = rnd.Rannyu();
		q2 = rnd.Rannyu();
		if(q1 < probability[2])	{
			this->Shift(ind1,rnd,cities);
			m_popolaz_new[contatore] = m_popolaz[ind1];
			contatore += 1;
			if(contatore == m_dim)	break;
		}
		if(q2 < probability[2])	{
			this->Shift(ind2,rnd,cities);
			m_popolaz_new[contatore] = m_popolaz[ind2];
			contatore += 1;
			if(contatore == m_dim)	break;
		}

Generaz<<"			Inversion"<<endl;
		q1 = rnd.Rannyu();
		q2 = rnd.Rannyu();
		if(q1 < probability[3])	{
			this->Inversion(ind1,rnd,cities);
			m_popolaz_new[contatore] = m_popolaz[ind1];
			contatore += 1;
			if(contatore == m_dim)	break;
		}

		if(q2 < probability[3])	{
			this->Inversion(ind2,rnd,cities);
			m_popolaz_new[contatore]= m_popolaz[ind2];
			contatore += 1;
			if(contatore == m_dim)	break;
		}

	 }
	// 3.2) Copio m_popolaz_new su m_popolaz
	this->Copy_popolaz();

	//4) Sorting
	this->sorting();
}

// fa la copia di m_popolaz_new su m_popolaz
void Popolazione :: Copy_popolaz(){
	//this->popolaz_erase();
	for( int i = 0 ; i < m_dim; i++){
		m_popolaz[i] = m_popolaz_new[i]; 			// basta cosi??
	}
	return;
}

// cancella tutto da m_populaz
void Popolazione ::  popolaz_erase(){
cout<<"erase 1"<<endl;
	m_popolaz.clear();
cout<<"erase 2"<<endl;

	return;
}     

// cancella tutto da m_populaz_copy
void Popolazione :: popolaz_new_erase(){
cout<<"erase 1"<<endl;
	m_popolaz_new.clear();
cout<<"erase 2"<<endl;
	return;
} 

