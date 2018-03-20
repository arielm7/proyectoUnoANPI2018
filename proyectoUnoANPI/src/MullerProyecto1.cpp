
#include <iostream>
#include <complex>
#include "math.h"
#include <boost/math/tools/polynomial.hpp>
using namespace boost::math::tools;
using namespace boost::math;

template <typename T>
T mullerMetod (const polynomial<T>& poly,T xr){
	float h = 1;
	float eps = 0.0001;
	int maxit = 10000;
	T x2 = xr;
	T x1 = xr + h;
	T x0 = xr - h;
	int iter = 0;
	T h0 = 0;
	T h1 = 0;
	T d0 = 0;
	T d1 = 0;
	T a = 0;
	T b = 0;
	T c = 0;
	T rad = 0;
	T den = 0;
	T dxr = 0;
	h0 = x1 - x0;
	h1 = x2 - x1;
	d0 = (poly.evaluate(x1) - poly.evaluate(x0))/h0;
	d1 = (poly.evaluate(x2) - poly.evaluate(x1))/h1;
	a = (d1 - d0)/ (h1 + h0);
	b = a*h1 + d1;
	c = poly.evaluate(x2);
	std::cout<<rad<< std::endl;
	if (isnan(rad)){
		std::cout << "Tenemos una raíz compleja" << std::endl;
		return NAN;
	}
	if (fabs(b+rad) > fabs(b-rad)){
		den = b+rad;		
	}
	if (fabs(b+rad) < fabs(b-rad)){
		den = b-rad;		
	}
	dxr = -2*c / den;
	xr = x2 + dxr;
	std::cout<<"La raiz encontrada es: " << xr << std::endl;
	return xr;
}

int main(){
	std::cout << "hola" << std::endl;
	double point = 1.0;
	polynomial<double> const poly{{-18,9,7,1,1}};
	double raiz = mullerMetod(poly, point);
	//std::cout << raiz << std::endl;
}
