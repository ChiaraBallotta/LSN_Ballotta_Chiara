#ifndef __INPUT__
#define __INPUT__

/////////////////////// VARIABILI UTILIZZATE //////////////////////
const double probability[5]={
    0.1,    //Pair_permutaz
    0.2,    //permutaz_N
    0.1,    //Shift
    0.2,    //Inversion
    0.6,    //CrossOver
};

const int N_cities = 34;      // numero di città contenute in ciascun individuo
const int N_individui = 3000;  // numero di individui contenuti in popolazione
const int N_generaz = 5000;   // numero di cicli in cui è aggiornata la popolaz mediante mutate

#endif