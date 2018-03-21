/*
 * laguerre.h
 *
 *  Created on: 12 mar. 2018
 *      Author: acacia
 */

#ifndef METODOS_LAGUERRE_LAGUERRE_H_
#define METODOS_LAGUERRE_LAGUERRE_H_

#include "../../Complejos/Complejo.h"
#include <boost/math/tools/polynomial.hpp>
#include <iostream>
#include <limits>

#define EPSS 1.0e-7
#define MR 8
#define MT 10
#define MAXIT (MT*MR)



//Hay que tener ciertas en cosas a la hora de utilizar el metodo de laguerre para el calculo de raices de ecuanciones.
//1. El metodo de Laguerre puede contener raices tanto reales, como complejas.
//2. Cuando los coeficientes de los polinomios son reales, y las raices complejas son las respuestas
//Estas aparecen por pares tanto el complejo como su conjugado.
//3. En cambio cuando el coeficiente es complejo, la raiz compleja no esta necesariamente relacionada.


template<class T>
class laguerre {
private:
	boost::math::tools::polynomial<T>* polinomio;
	int degree;
	T* x;
	int* iteraciones;


public:
	laguerre();
	void initLaguerre(boost::math::tools::polynomial<T>* polinomio, double* x);
	void initLaguerre(boost::math::tools::polynomial<T>* polinomio, float* x);
	void initLaguerre(boost::math::tools::polynomial<T>* polinomio, std::complex<float>* x);
	void initLaguerre(boost::math::tools::polynomial<T>* polinomio, std::complex<double>* x);
	void laguer(boost::math::tools::polynomial<std::complex<float>> polinomio, int grado,std::complex<float> x);
	void laguer(boost::math::tools::polynomial<std::complex<double>> polinomio,int grado, std::complex<double> x);
	void convPolinomioComplejo(boost::math::tools::polynomial<float> polinomio, float* x);
	void convPolinomioComplejo(boost::math::tools::polynomial<double> polinomio, double* x);
	boost::math::tools::polynomial<complex<float>> zroots(boost::math::tools::polynomial<complex<float>> polinomio, bool polish);
	boost::math::tools::polynomial<complex<double>> zroots(boost::math::tools::polynomial<complex<double>> polinomio, bool polish);
	virtual ~laguerre();

/**
 * Aqui se encuentran los get de la clase.
 */
public:
	boost::math::tools::polynomial<T>* getPolinomio();
	int getDegree();
	T* getX();
	void setPolinomio();
	void setDegree();
	void setX();
};

/**
 * Constructor de la clase se encarga de construir la clase
 */
template<class T>
laguerre<T>::laguerre(){
	this->degree=0;
	this->polinomio = NULL;
	this->x = NULL;
	this->iteraciones=0;
}

/**
 * Destructor de la clase
 */
template<class T>
laguerre<T>::~laguerre() {
	// TODO Auto-generated destructor stub
}


template<class T>
void laguerre<T>::initLaguerre(boost::math::tools::polynomial<T>* polinomio, double* x){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	convPolinomioComplejo(*polinomio,x);
}

template<class T>
void laguerre<T>::initLaguerre(boost::math::tools::polynomial<T>* polinomio, float* x){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	convPolinomioComplejo(*polinomio,x);
}

template<class T>
void laguerre<T>::initLaguerre(boost::math::tools::polynomial<T>* polinomio, std::complex<double>* x){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
}

template<class T>
void laguerre<T>::initLaguerre(boost::math::tools::polynomial<T>* polinomio, std::complex<float>* x){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
}

template<class T>
void laguerre<T>::laguer(boost::math::tools::polynomial<std::complex<float>> polinomio, int grado ,std::complex<float> x){
	Complejo<float>* AC = new Complejo<float>();
	int iter, j;
	float abx,abp,abm,err;
	std::complex<float> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	for(iter = 1; iter<=MAXIT;iter++){
		*this->iteraciones = iter;
		b=polinomio[grado];
		d=f= std::complex<T>(0.0,0.0);
		abx = AC->abs(x);

		for(j=grado-1;j>=0;j--){
			f=((x*f)+d);
			d=((x*d)+b);
			b=((x*b)+polinomio[j]);
			err=AC->abs(b)*abx*err;
		}
		err*= EPSS;

		if(AC->abs(b) <= err) return;
		g= d/b;
		g2 = g*g;
		h=(g2-(T(2)*(f/b)));
		sq = AC->sqrt(((float)(grado-1))*(((float)(grado)*h)-g2));
		gp=g+sq;
		gm=g-sq;
		abp=AC->abs(gp);
		abm=AC->abs(gm);
		if(abp<abm) gp=gm;
		dx=((fmax(abp,abm) > 0.0 ? ((std::complex<float>(grado,0.0))/gp) : (1+abx)*(std::complex<float>(cos(iter),sin(iter)))));
		x1=x-dx;
		if(x == x1) return ;
		if(iter%MT) x=x1;
		else{
			x = x - AC->cmul(frac[iter/MT],dx);
		}

	}

	std::cout << "Ocurrio un error" <<std::endl;
	return;
}

template<class T>
void laguerre<T>::laguer(boost::math::tools::polynomial<std::complex<double>> polinomio, int grado ,std::complex<double> x){
	Complejo<double>* AC = new Complejo<double>();
	int iter, j;
	double abx,abp,abm,err;
	std::complex<double> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	for(iter = 1; iter<=MAXIT;iter++){
		*this->iteraciones = iter;
		b=polinomio[grado];
		d=f= std::complex<T>(0.0,0.0);
		abx = AC->abs(x);//sqrt(pow(this->x->real(),2) + pow(this->x->imag(),2));

		for(j=grado-1;j>=0;j--){
			f=((x*f)+d);
			d=((x*d)+b);
			b=((x*b)+polinomio[j]);
			err=AC->abs(b)*abx*err;
		}
		err*= EPSS;

		if(AC->abs(b) <= err) return;
		g= d/b;
		g2 = g*g;
		h=(g2-(T(2)*(f/b)));
		sq = AC->sqrt(((double)(grado-1))*(((double)(grado)*h)-g2));
		gp=g+sq;
		gm=g-sq;
		abp=AC->abs(gp);
		abm=AC->abs(gm);
		if(abp<abm) gp=gm;
		dx=((fmax(abp,abm) > 0.0 ? ((std::complex<double>(grado,0.0))/gp) : (1+abx)*(std::complex<double>(cos(iter),sin(iter)))));
		x1=x-dx;
		if(x == x1) return ;
		if(iter%MT) x=x1;
		else x= x - AC->cmul(frac[iter/MT],dx);

	}

	std::cout << "La cantidad de iteraciones son muchas" <<std::endl;
	return;
}

template<class T>
void laguerre<T>::convPolinomioComplejo(boost::math::tools::polynomial<double> polinomio, double* x){
	boost::math::tools::polynomial<std::complex<double>>  elementos;
	for(int i = 0; i<= this->degree; i++){
		elementos[i] = complex<double>(polinomio[i],0.0);
	}

}

template<class T>
void laguerre<T>::convPolinomioComplejo(boost::math::tools::polynomial<float> polinomio, float* x){
	boost::math::tools::polynomial<std::complex<float>> elementos;
	for(int i = 0; i<= this->degree; i++){
			elementos[i] = complex<float>(polinomio[i],0.0);
	}
}

template<class T>
boost::math::tools::polynomial<std::complex<double>> laguerre<T>::zroots(boost::math::tools::polynomial<std::complex<double>> polinomio ,bool polish){
	Complejo<double> * AC = new Complejo<double>();
	int i, j,jj;
	complex<double> x,b,c;
	complex<double> roots[polinomio.size()];
	for (i=0;i<polinomio.size();i++){
        roots[i]=new complex<double>(0.0,0.0);
	}
	complex<double> ad[MAXIT];

	for(j=0; j<MAXIT; j++){
		ad[j]=new complex<double>(0.0,0.0);
	}

	for (j=0;j<=polinomio.degree();j++) ad[j]=polinomio[j];

	for (j=polinomio.degree();j>=1;j--) {
		x = new complex<double>(0.0,0.0);
		laguer(ad,j, x);
		if (abs(x.imag()) <= 2.0*numeric_limits<T>::epsilon()*abs(x.real())) x.imag()=0.0;
		roots[j]=x;
		b=ad[j];
		for (jj=j-1;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=AC->add(AC->mul(x,b),c);
		}
	}

	if (polish)
		for (j=1;j<=polinomio.degree();j++)
			laguer(polinomio,polinomio.degree(),roots[j]);
	for (j=2;j<=polinomio.degree();j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].real() <= x.real()) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
  return roots;
}

template<class T>
boost::math::tools::polynomial<std::complex<float>> laguerre<T>::zroots(boost::math::tools::polynomial<std::complex<float>> polinomio ,bool polish){
	Complejo<float> * AC = new Complejo<float>();
	int i, j,jj;
	complex<float> x,b,c;
	complex<float> roots[polinomio.size()];
	for (i=0;i<polinomio.size();i++){
        roots[i]=new complex<float>(0.0,0.0);
	}
	complex<double> ad[MAXIT];

	for(j=0; j<MAXIT; j++){
		ad[j]=new complex<float>(0.0,0.0);
	}

	for (j=0;j<=polinomio.degree();j++) ad[j]=polinomio[j];

	for (j=polinomio.degree();j>=1;j--) {
		x = new complex<float>(0.0,0.0);
		laguer(ad,j, x);
		if (abs(x.imag()) <= 2.0*numeric_limits<T>::epsilon()*abs(x.real())) x.imag()=0.0;
		roots[j]=x;
		b=ad[j];
		for (jj=j-1;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=AC->add(AC->mul(x,b),c);
		}
	}

	if (polish)
		for (j=1;j<=polinomio.degree();j++)
			laguer(polinomio,polinomio.degree(),roots[j]);
	for (j=2;j<=polinomio.degree();j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].real() <= x.real()) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
  return roots;
}

#endif /* METODOS_LAGUERRE_LAGUERRE_H_ */
