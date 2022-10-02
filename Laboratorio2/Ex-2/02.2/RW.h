
#ifndef __RW__
#define __RW__

#include "random.h"
/////////// RW DISCRETO ////////
class Posizione_D{
	public:
		Posizione_D(Random&);
		Posizione_D(int , int , int , Random& );
        ~Posizione_D();
		//metodi
		int GetX();
		int GetY();		
		int GetZ();
		int Sign(Random&);
		int Direction(Random&);
		double Dist2();
		void Update(int[],Random&);		
	private: 
		int m_x,m_y,m_z;
};

/////////// RW CONTINUA ////////
class Posizione_C{
	public:
		Posizione_C(Random&,int);
		Posizione_C(double , double , double , Random&,int );
        ~Posizione_C();
		//metodi
		double GetX();
		double GetY();		
		double GetZ();
		double Dist2();
		void Update(double[],Random&);		
	private: 
		double m_x,m_y,m_z;
		int m_a;		// lunghezza dello step
};



#endif