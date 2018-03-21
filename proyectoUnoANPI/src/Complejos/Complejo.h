/*
 * Complejo.h
 *
 *  Created on: 20 mar. 2018
 *      Author: acacia
 */
#ifndef COMPLEJOS_COMPLEJO_H_
#define COMPLEJOS_COMPLEJO_H_

#include <math.h>
#include <complex.h>
#include <complex>

using namespace std;

template<class T>
class Complejo : public complex<T>{
public:
	complex<T> add(complex<T> a , complex<T> b);
	complex<T> sub(complex<T> a , complex<T> b);
	complex<T> mul(complex<T> a , complex<T> b);
	complex<T> conj(complex<T>  z);
	complex<T> div(complex<T> a , complex<T> b);
	T abs(complex<T> z);
	complex<T> sqrt(complex<T> z);
	complex<T> cmul(T data, complex<T> z);
	Complejo();
	virtual ~Complejo();
};

template<class T>
complex<T> Complejo<T>::add(complex<T> a, complex<T> b){
	return a + b;
}

template<class T>
complex<T> Complejo<T>::sub(complex<T> a, complex<T> b){
	return b-a;
}

template<class T>
complex<T> Complejo<T>::mul(complex<T> a, complex<T> b){
	return a*b;
}

template<class T>
complex<T> Complejo<T>::conj(complex<T> z){
	return new complex<T>(z.real(),-z.imag());
}

template<class T>
complex<T> Complejo<T>::div(complex<T> a, complex<T> b){
	return a/b;
}

template<class T>
T Complejo<T>::abs(complex<T> z){
	T x,y,ans,temp;
	x=abs(z.real());
	y=abs(z.imag());
	if (x == 0.0)
		ans=y;
	else if (y == 0.0)
		ans=x;
	else if (x > y) {
		temp=y/x;
		ans=x*sqrt(1.0+temp*temp);
	} else {
		temp=x/y;
		ans=y*sqrt(1.0+temp*temp);
	}
	return ans;
}


template<class T>
complex<T> Complejo<T>::cmul(T data, complex<T> z){
	return new complex<T>(data*z.real(), data*z.imag());
}

template<class T>
complex<T> Complejo<T>::sqrt(complex<T> z){
	T real, imag;
	T x,y,w,r;
	if ((z.real() == 0.0) && (z.imag() == 0.0)) {
		return new complex<T>();
	} else {
		x=abs(z.real());
		y=abs(z.imag());
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r >= 0.0) {
			real=w;
			imag=z.imag()/(2.0*w);
		} else {
			imag=(z.imag() >= 0) ? w : -w;
			real=z.imag()/(2.0*imag);
		}
		return new complex<T>(real, imag);
	}
}

#endif /* COMPLEJOS_COMPLEJO_H_ */
