
#include<vector>
#include <iostream>
#include "random.h"
#include "input.h"


using namespace std;


#ifndef __INDIVIDUO__
#define __INDIVIDUO__

//////////////  CLASSE INDIVIDUO  /////////////////////////////////
class Individuo {

public:
//generatore Random
  Individuo();
  Individuo(int ,Random&);
  Individuo(int*,int);
  ~Individuo();

  void Individuo_new(Random&);

  //vector<int> m_ind;	//vettore che contiene la sequenza

//double length(mat);
  void Set_length(double**);
  double Get_length();
  void Modif_city(int,int);
  int Get_cities(int);
  int Get_n();
  int Get_size();


private:
  double m_length;
  vector<int> m_ind;	//vettore che contiene la sequenza
  int m_n;            // Dimensione vector Individuo--> numero di citta visitate

};

#endif // __INDIVIDUO__

#ifndef __POPOLAZIONE__
#define __POPOLAZIONE__
/////////////////////   CLASSE POPOLAZIONE    ////////////////////////////

class Popolazione {

public:
  Popolazione(int, int, Random& ,double**);
  ~Popolazione();
  int Get_dim();
  void sorting();
  int selez(Random&);
  Individuo Get_individuo(int);       //restituisce l'individuo i-esimo dal bach della popolazione
  void Set_individuo(Individuo,int);    


  //Permutazioni
  void Pair_permutaz(int,Random&,double**);
  void permutaz_N(int,Random&,double**);
  void Shift(int,Random&,double**);
  void Inversion(int,Random&,double**);
  void Crossover(int, int,Random&,double**);

  // Mutazione
  const double *Prob = probability; // basic mutation probabilit
  void Mutation(Random&, double**,ofstream&);

  void Copy_popolaz();        // fa la copia di m_popolaz su m_popolaz_copy
  void popolaz_erase();       // cancella tutto da m_populaz
  void popolaz_new_erase();  // cancella tutto da m_populaz_copy  

private:
  int m_dim;    //numero di individui contenuti all' interno della popolazione
  vector<Individuo> m_popolaz;  //vettore che conterrà la popolazione
  vector<Individuo> m_popolaz_new;  // vettore che contiene la nuova popolazione dopo mutazione
  
};

#endif // __POPOLAZIONE__

////////////// FUNZIONI /////////////////////////
double dist(double,double,double,double);
bool check(Individuo);
void scambiaByRef(Individuo*,Individuo*);
