#ifndef __funzioniR1_h__
#define __funzioniR1_h__

#include <iostream>
#include <fstream>
#include <string>
#include<cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include "math.h"

using namespace std;

////////////Funzioni base in R1:

		//Funzione Base
		
		//COMBINAZ3FUNZ : Combinazione lineare di tre funzioni base in R1

		//COSENO

		//PARABOLA

		//RETTA

		//SENO

		//GAUSSIANA

		//LOGARITMO

		//x^(a)

		//DIVISIONE


///////////////////////////FUNZIONE BASE

class FunzioneBase{
  public:
  	virtual double Eval(double x) const = 0;
};

///////////////////////////Permette il calcolo della combinazione Lineare di 3 funzioni: CombinazLineare = d*a(x) + e*b(x) + f*c(x) + g*a(x)*b(x) + h*a(x)*c(x) + i*c(x)*b(x)

class Combinaz3Funz :public FunzioneBase{

  public:
  
  	Combinaz3Funz(FunzioneBase * a, FunzioneBase * b,FunzioneBase* c,double d, double e ,double f,double g, double h, double i){
	   m_a = a; m_b = b;m_c = c;m_d = d; m_e = e;m_f = f;m_g = g;m_h = h;m_i = i;
	  }

  	~Combinaz3Funz(){};

  	virtual double Eval(double x) const{

	  	return m_d*(m_a->Eval(x))+m_e*(m_b->Eval(x))+m_f*(m_c->Eval(x))+m_g*(m_a->Eval(x))*(m_b->Eval(x))+ m_h*(m_a->Eval(x))*(m_c->Eval(x))+m_i*(m_c->Eval(x))*(m_b->Eval(x));
 	}  

  private:
  	double m_d, m_e,m_f,m_g,m_h,m_i;
  	FunzioneBase * m_a;
  	FunzioneBase * m_b;
 	FunzioneBase * m_c;
};

///////////////////////////Funzione COSENO y = a * cos(w*x + c)

class Coseno:public FunzioneBase{
 public:
 
	Coseno(double a, double w,double c){
	  m_a = a; m_w = w;m_c= c;
	  }

	~Coseno(){};

	virtual double Eval(double x)const{
	  	return m_a*cos(m_w*x + m_c);
	  }

//creo l' overloading dell' operatore ()
	double operator()(double x)const{return m_a*cos(m_w*x + m_c);};

	void SetA(double a) {m_a = a;}
	void SetW(double w) {m_w = w;}
	void SetC(double c) {m_c = c;}

	double GetA() const {return m_a;}
	double GetW() const {return m_w;}
	double GetC() const {return m_c;}

//implementiamo la copia di elemento
	Coseno(const Coseno& V){
		m_a = V.GetA();
		m_w = V.GetW();
		m_c = V.GetC();};
	Coseno& operator=(const Coseno& V){
		m_a = V.GetA();
		m_w = V.GetW();
		m_c = V.GetC();
		return *this;
	};
private:
  double m_a, m_w,m_c;
};

/////////////////////////// PARABOLA y = ax^2 + bx + c

class Parabola:public FunzioneBase{
  public:
  	Parabola(){
		m_a = 0;
	  	m_b = 0;
	  	m_c = 0;
	}
  	Parabola(double a, double b, double c){
		m_a = a;
	  	m_b = b;
	  	m_c = c;
	}
  	~Parabola(){};
  	virtual double Eval(double x) const{return m_a*x*x + m_b*x + m_c;}

  //creo l' overloading dell' operatore ()
  	double operator()(double x)const{return m_a*x*x + m_b*x+m_c;}
  	void SetA(double a) {m_a = a;}
  	void SetB(double b) {m_b = b;}
  	void SetC(double c) {m_c = c;}

  	double GetA() const {return m_a;}
  	double GetB() const {return m_b;}
  	double GetC() const {return m_c;}

  //implementiamo la copia di elemento
  	Parabola(const Parabola& V){
		m_a = V.GetA();
		m_b = V.GetB();
		m_c = V.GetC();
	}
  	Parabola& operator=(const Parabola& V){
		m_a = V.GetA();
		m_b = V.GetB();
		m_c = V.GetC();
		return*this;
	}

private:
  	double m_a, m_b, m_c;
};

/////////////////////////// RETTA y = a* x + b

class Retta:public FunzioneBase{
  public:

  	Retta(double a, double b){m_a = a; m_b = b;}
  	~Retta(){};

  	virtual double Eval(double x) const{return m_a*x + m_b;}

  //creo l' overloading dell' operatore ()
  	double operator()(double x)const{return m_a*x + m_b;};

  	void SetA(double a) {m_a = a;}
  	void SetB(double b) {m_b = b;}

  	double GetA() const {return m_a;}
  	double GetB() const {return m_b;}
  
  //implementiamo la copia di elemento
	Retta(const Retta& V){
		m_a = V.GetA();
		m_b = V.GetB();
	};
	Retta& operator=(const Retta& V){
		m_a = V.GetA();
		m_b = V.GetB();	
		return *this;
	};
private:
  double m_a, m_b;
};

///////////////////////////Funzione SENO y = a * sin(w*x + c)

class Seno:public FunzioneBase{
  public:
  
  	Seno(double a, double w,double c){m_a = a; m_w = w;m_c= c;}
  	~Seno(){};

  	virtual double Eval(double x) const{return m_a*sin(m_w*x + m_c);}
  //creo l' overloading dell' operatore ()
  	double operator()(double x)const{return m_a*sin(m_w*x + m_c);}

 	void SetA(double a) {m_a = a;}
  	void SetW(double w) {m_w = w;}
  	void SetC(double c) {m_c = c;}

  	double GetA() const {return m_a;}
  	double GetW() const {return m_w;}
  	double GetC() const {return m_c;}

  //implementiamo la copia di elemento

	Seno(const Seno& V){
		m_a = V.GetA();
		m_w = V.GetW();
		m_c = V.GetC();};
	Seno& operator=(const Seno& V){
		m_a = V.GetA();
		m_w = V.GetW();
		m_c = V.GetC();
		return *this;
	};
private:
  double m_a, m_w,m_c;
};


/////////////////////////////////////////////////////////////////Gaussiana
class Gaussiana:public FunzioneBase{
public:

  Gaussiana(double sigma, double mu){m_sigma = sigma; m_mu = mu;}
  ~Gaussiana(){};
  virtual double Eval(double x) const{
	  double exp = (-1)*(pow(x-m_mu,2))/(2*m_sigma*m_sigma);
	  return (1/(m_sigma*sqrt(2*M_PI)))*pow(M_E,exp);
	  }

  void SetS(double sigma) {m_sigma = sigma;}
  void SetM(double mu) {m_mu = mu;}
  double GetS() const {return m_sigma;}
  double GetM() const {return m_mu;}
  double GetMax() const{return Eval(m_mu);}
  


private:
  double m_sigma, m_mu;
};

///////////////////////////Funzione LOGARITMO y = a * log(b)

class Logaritmo:public FunzioneBase{
  public:
  
  	Logaritmo(double a, FunzioneBase* b){m_a = a; m_b = b;}
  	~Logaritmo(){};

  	virtual double Eval(double x) const{return m_a*log(m_b->Eval(x));}
 
private:
  double m_a;
  FunzioneBase* m_b;
};
///////////////////////////Funzione POTENZA y = a + d*( x + b)^(c)

class Potenza:public FunzioneBase{
  public:
  
  	Potenza(double a,double b,double c,double d,FunzioneBase* x){m_a = a; m_b = b;m_c = c;m_x = x;m_d = d;}
  	~Potenza(){};

  	virtual double Eval(double x) const{return m_a + m_d * pow(m_x->Eval(x) + m_b, m_c);}
 
private:
  double m_a, m_b, m_c,m_d;
  FunzioneBase* m_x;
};

///////////////////////////Funzione DIVISIONE y = a/b*x

class Divisione:public FunzioneBase{
  public:
  
  	Divisione(double a,double b,FunzioneBase* x){m_a = a; m_b = b;m_x = x;}
  	~Divisione(){};

  	virtual double Eval(double x) const{
		  if(m_x->Eval(x) == 0) {
			  cout<<"Errore, il divisore si annulla in x = "<<x<<endl;
			  exit(-1);
		  }
		  return m_a/(m_b * m_x->Eval(x));}
 
private:
  double m_a, m_b;
  FunzioneBase* m_x;
};
#endif
