/*
 * laguerre.h
 *
 *  Created on: 12 mar. 2018
 *      Author: acacia
 */

#ifndef METODOS_LAGUERRE_DEFLATEDLAGUERRE_H_
#define METODOS_LAGUERRE_DEFLATEDLAGUERRE_H_


#include "../../deflateAlgorithms.h"
#include <iostream>
#include <boost/math/tools/polynomial.hpp>
#include <complex>
#include <complex.h>
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

using namespace std;

template<class T>
class deflatedLaguerre {
private:
	boost::math::tools::polynomial<T>* polinomio;
	int degree;
	T x;
	int iteraciones;
	bool pulido;


public:
	deflatedLaguerre();
	void initLaguerre(boost::math::tools::polynomial<double>* polinomio, double &x, bool pulido);
	void initLaguerre(boost::math::tools::polynomial<float>* polinomio, float &x, bool pulido);
	void initLaguerre(boost::math::tools::polynomial<std::complex<float>>* polinomio, std::complex<float> &x, bool pulido);
	void initLaguerre(boost::math::tools::polynomial<std::complex<double>>* polinomio, std::complex<double> &x, bool pulido);
	void laguer(boost::math::tools::polynomial<std::complex<float>> polinomio, int grado,std::complex<float> &x);
	void laguer(boost::math::tools::polynomial<std::complex<double>> polinomio,int grado, std::complex<double> &x);
	void convPolinomioComplejo(boost::math::tools::polynomial<float> polinomio, complex<float> &x);
	void convPolinomioComplejo(boost::math::tools::polynomial<double> polinomio, complex<double> &x);
	vector<complex<float>> zroots(boost::math::tools::polynomial<complex<float>> polinomio, bool polish);
	vector<complex<double>> zroots(boost::math::tools::polynomial<complex<double>> polinomio, bool polish);
	virtual ~deflatedLaguerre();

/**
 * Aqui se encuentran los get de la clase.
 */
public:
	void setPulido(bool pulido);
	bool getPulido();

};

/**
 * Constructor de la clase se encarga de construir la clase
 */
template<class T>
deflatedLaguerre<T>::deflatedLaguerre(){
	this->degree=0;
	this->polinomio = NULL;
	this->x = T(0);
	this->iteraciones=0;
	this->pulido = false;
}

/**
 * Destructor de la clase
 */
template<class T>
deflatedLaguerre<T>::~deflatedLaguerre() {
}


template<class T>
void deflatedLaguerre<T>::initLaguerre(boost::math::tools::polynomial<double>* polinomio, double &x, bool pulido){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	this->pulido = pulido;
	complex<double> p(x,T(0));
	convPolinomioComplejo(*polinomio,p);
}

template<class T>
void deflatedLaguerre<T>::initLaguerre(boost::math::tools::polynomial<float>* polinomio, float &x, bool pulido){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	std::complex<float> nuevo(x,T(0));
	convPolinomioComplejo(*polinomio,nuevo);
}

template<class T>
void deflatedLaguerre<T>::initLaguerre(boost::math::tools::polynomial<std::complex<double>>* polinomio, std::complex<double> &x, bool pulido){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	zroots(*polinomio,this->pulido);
}

template<class T>
void deflatedLaguerre<T>::initLaguerre(boost::math::tools::polynomial<std::complex<float>>* polinomio, std::complex<float> &x, bool pulido){
	this->polinomio = polinomio;
	this->degree = this->polinomio->degree();
	this->x = x;
	zroots(*polinomio,this->pulido);
}

template<class T>
void deflatedLaguerre<T>::laguer(boost::math::tools::polynomial<std::complex<double>> polinomio, int grado ,std::complex<double> &x){
	int iter, j;
	double abx,abp,abm,err;
	std::complex<double> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static double frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
	for(iter = 1; iter<=MAXIT;iter++){
		this->iteraciones = iter;
		b=polinomio[grado];
		d=f= std::complex<double>(0.0,0.0);
		abx = abs(x);

		for(j=grado-1;j>=0;j--){
			f=((x*f)+d);
			d=((x*d)+b);
			b=((x*b)+polinomio[j]);
			err=abs(b)*abx*err;
		}
		err*= EPSS;

		if(abs(b) <= err) return;
		g= d/b;
		g2 = g*g;
		h=(g2-(2.0*(f/b)));
		sq = sqrt(((double)(grado-1))*(((double)(grado)*h)-g2));
		gp=g+sq;
		gm=g-sq;
		abp=abs(gp);
		abm=abs(gm);
		if(abp<abm) gp=gm;
		dx=((fmax(abp,abm) > 0.0 ? ((std::complex<double>(grado,0.0))/gp) : (1+abx)*(std::complex<double>(cos(iter),sin(iter)))));
		x1=x-dx;
		if(x == x1) return ;
		if(iter%MT) x=x1;
		else{
			x = x - (frac[iter/MT] * dx);
		}

	}

	std::cout << "No encontro la raiz" <<std::endl;
	return;
}

template<class T>
void deflatedLaguerre<T>::laguer(boost::math::tools::polynomial<std::complex<float>> polinomio, int grado ,std::complex<float> &x){
	int iter, j;
	float abx,abp,abm,err;
	std::complex<float> dx,x1,b,d,f,g,h,sq,gp,gm,g2;
	static float frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

	for(iter = 1; iter<=MAXIT;iter++){
		this->iteraciones = iter;
		b=polinomio[grado];
		d=f= std::complex<float>(0.0,0.0);
		abx = abs(x);

		for(j=grado-1;j>=0;j--){
			f=((x*f)+d);
			d=((x*d)+b);
			b=((x*b)+polinomio[j]);
			err=abs(b)*abx*err;
		}
		err*= EPSS;

		if(abs(b) <= err) return;
		g= d/b;
		g2 = g*g;
		h = (g2-((float)2 * (f/b)));
		sq = sqrt(((float)(grado-1))*(((float)(grado)*h)-g2));
		gp=g+sq;
		gm=g-sq;
		abp=abs(gp);
		abm=abs(gm);
		if(abp<abm) gp=gm;
		dx=((fmax(abp,abm) > 0.0 ? ((std::complex<float>(grado,0.0))/gp) : (1+abx)*(std::complex<float>(cos(iter),sin(iter)))));
		x1=x-dx;
		if(x == x1) return ;
		if(iter%MT) x=x1;
		else{
			x = x - (frac[iter/MT] * dx);
		}

	}

	std::cout << "No encontro la raiz" <<std::endl;
	return;
}

template<class T>
void deflatedLaguerre<T>::convPolinomioComplejo(boost::math::tools::polynomial<double> polinomio, complex<double> &x){
	std::complex<double> arreglo[polinomio.size()];
	for(int i = 0; i< (int)polinomio.size(); i++){
		arreglo[i] = complex<double>(polinomio[i],T(0));
	}
	boost::math::tools::polynomial<std::complex<double>>  elementos(arreglo, polinomio.degree());
	zroots(elementos,this->pulido);

}

template<class T>
void deflatedLaguerre<T>::convPolinomioComplejo(boost::math::tools::polynomial<float> polinomio, complex<float> &x){
	std::complex<float> arreglo[polinomio.size()];
		for(int i = 0; i< (int)polinomio.size(); i++){
			arreglo[i] = complex<float>(polinomio[i],T(0));
		}
		boost::math::tools::polynomial<std::complex<float>>  elementos(arreglo, polinomio.degree());
		zroots(elementos,this->pulido);
}

template<class T>
vector<complex<double>> deflatedLaguerre<T>::zroots(boost::math::tools::polynomial<std::complex<double>> polinomio ,bool polish){
	int i, j;
	complex<double> x,b,c, residuo;
	std::vector<complex<double>> roots(polinomio.degree());
	boost::math::tools::polynomial<std::complex<double>> ad = polinomio;
	boost::math::tools::polynomial<std::complex<double>> residuo2({0,0});
	for (j=polinomio.degree();j>=1;j--) {
		x = complex<double>(0.0,0.0);
		residuo = complex<double>(0.0,0.0);
		laguer(ad,j, x);
		if (abs(x.imag()) <= (double)2*numeric_limits<double>::epsilon()*abs(x.real())){
			x.imag(0);
			ad = deflate<std::complex<double>>(ad,x,residuo);

		}else{
			if(abs(x.real()) <= (double)2*numeric_limits<double>::epsilon()*abs(x.imag())){
				x.real(0);
			}
			if(typeid(T) == typeid(double)){
				cout << "Se encontro una raiz imaginaria" << endl;
			}
			ad=deflate2Laguerre<std::complex<double>>(ad,x,residuo2);
			roots[j-1]=std::conj(x);

			j--;
		}
		roots[j-1] = x;

	}

	//Proceso de pulido de raices.
	if (polish)
		for (j=1;j<=(int)polinomio.degree();j++)
			laguer(polinomio,(int)polinomio.degree(),roots[j]);
	for (j=2;j<=(int)polinomio.degree();j++) {
		x=roots[j];
		for (i=j-1;i>=1;i--) {
			if (roots[i].real() <= x.real()) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}

	std::cout<<"raices encontradas: "<<std::endl;
	for(int i=0;i<(int)polinomio.degree();i++ ){
		std::cout<<roots[i]<<std::endl;
	}
	return roots;
}



template<class T>
vector<complex<float>> deflatedLaguerre<T>::zroots(boost::math::tools::polynomial<std::complex<float>> polinomio ,bool polish){
	int i, j;
	complex<float> x,b,c, residuo;
	std::vector<complex<float>> roots(polinomio.degree());
	boost::math::tools::polynomial<std::complex<float>> ad = polinomio;
	boost::math::tools::polynomial<std::complex<float>> residuo2({0,0});
	for (j=polinomio.degree();j>=1;j--) {
		x = complex<float>(0.0,0.0);
		residuo = complex<float>(0.0,0.0);
		laguer(ad,j, x);
		if (abs(x.imag()) <= (float)2*numeric_limits<float>::epsilon()*abs(x.real())){
			x.imag(0);
			ad = deflate<std::complex<float>>(ad,x,residuo);
		}else{
			if(abs(x.real()) <= (float)2*numeric_limits<float>::epsilon()*abs(x.imag())){
				x.real(0);
			}
			if(typeid(T) == typeid(float)){
				cout << "Se encontro una raiz imaginaria" << endl;
			}
			ad=deflate2Laguerre<std::complex<float>>(ad,x,residuo2);
			roots[j-1]=std::conj(x);
			j--;
		}
		roots[j-1] = x;
	}
	//Proceso de pulido de raices.
	if (polish)
			for (j=1;j<=(int)polinomio.degree();j++)
				laguer(polinomio,(int)polinomio.degree(),roots[j]);
		for (j=2;j<=(int)polinomio.degree();j++) {
			x=roots[j];
			for (i=j-1;i>=1;i--) {
				if (roots[i].real() <= x.real()) break;
				roots[i+1]=roots[i];
			}
			roots[i+1]=x;
		}
	std::cout<<"raices encontradas: "<<std::endl;
	for( i=0;(int)i<(int)polinomio.degree();i++ ){
		std::cout<<roots[i]<<std::endl;
	}
	return roots;
}

#endif /* METODOS_LAGUERRE_DEFLATEDLAGUERRE_H_ */
