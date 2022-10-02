#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include "RW.h"

using namespace std;

/////////// RW DISCRETO ////////
Posizione_D ::Posizione_D(Random& rand){
	m_x = 0;
	m_y = 0;
	m_z = 0;
}
Posizione_D :: Posizione_D(int x, int y, int z, Random& rand){
	m_x = x;
	m_y = y;
	m_z = z;
}


int Posizione_D :: GetX(){
	return m_x;
}

int Posizione_D ::GetY(){
	return m_y;
}
		
int Posizione_D ::GetZ(){
	return m_z;
}
		
//Genera ka direzione in cui il Random-Walk discreto si muove
int Posizione_D :: Sign(Random& rand){
	double x;
	x = rand.Rannyu(0,1);
	if(x <= 0.5 )	return -1;
	if(x > 0.5) return 1;
	return 0;	// c'è stato un problema
}

//Genera la Direzione in cui il Random-Walk discreto si muove
int Posizione_D :: Direction(Random& rand){
    double x;
	x = rand.Rannyu(0,3);
	if(x <= 1 )	return 1;			//x
	if( 1 < x && x <= 2) return 2;	//y
	if( x > 2 && x <= 3) return 3;	//z
	return 0;	// c'è stato un problema

}

double Posizione_D :: Dist2(){
	int x = m_x;
	int y = m_y;
	int z = m_z;
	double dist2 =  x*x + y*y + z*z;
	return dist2;
}

//metodo che aggiorna la posizione del random walk
void Posizione_D :: Update(int pos[3],Random& rand){
	int dir = this->Direction(rand);
	int sign = this->Sign(rand);			
	if ( dir == 1)	pos[0] = pos[0] + sign;
	if ( dir == 2)	pos[1] = pos[1] + sign;
	if ( dir == 3)	pos[2] = pos[2] + sign;				
}

/////////// RW CONTINUO ////////
Posizione_C ::Posizione_C(Random& rand,int a){
	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_a = a;
}
Posizione_C :: Posizione_C(double x, double y, double z, Random& rand,int a){
	m_x = x;
	m_y = y;
	m_z = z;
	m_a = a;
}


double Posizione_C :: GetX(){
	return m_x;
}

double Posizione_C ::GetY(){
	return m_y;
}
		
double Posizione_C ::GetZ(){
	return m_z;
}
		

double Posizione_C :: Dist2(){
	double x = m_x;
	double y = m_y;
	double z = m_z;
	double dist2 =  x*x + y*y + z*z;
	return dist2;
}

//metodo che aggiorna la posizione del random walk
void Posizione_C :: Update(double pos[3],Random& rand){
	double theta = rand.Rannyu(0,M_PI);
	double phi = rand.Rannyu(0,2*M_PI);			
	pos[0] = pos[0] + m_a* sin(theta)*cos(phi);
	pos[1] = pos[1] + m_a*sin(theta)*sin(phi);
	pos[2] = pos[2] + m_a*cos(phi);			
}
