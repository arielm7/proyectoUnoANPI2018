/*
 * deflateAlgorithms.h
 *
 *  Created on: 21 mar. 2018
 *      Author: acacia
 */

#ifndef DEFLATEALGORITHMS_H_
#define DEFLATEALGORITHMS_H_

#include <boost/math/tools/polynomial.hpp>
#include <complex>
#include <iostream>


template<typename T>
boost::math::tools::polynomial<T> deflate(const boost::math::tools::polynomial<T>  poly, const T& root, T& residuo){

	boost::math::tools::polynomial<T>   result=poly;
	int n=poly.size()-1;
	residuo=poly[n];
	result[n]=T(0);
	for(int i=n-1;i>=0;i--) {

		T swap=result[i];
		result[i]=residuo;
		residuo=swap+residuo*root;


		}
	return result;
}//end deflate

template<typename T>
boost::math::tools::polynomial<T>  deflate2Laguerre(const boost::math::tools::polynomial<T> & poly, const
		T& root,boost::math::tools::polynomial<T> & residuo){



	T c=(root.real()*root.real())+(root.imag()*root.imag());

	T b= (T)((double)-2*root.real());
	T a=1;
	boost::math::tools::polynomial<T>  result=poly;
	boost::math::tools::polynomial<T>  cuadratic={{c,b,a}};
	int n=poly.size()-1;
	residuo=poly;

	int k,j;
	for (j=0;j<=n;j++) {

	result[j]=T(0);
	}

	for (k=n-2;k>=0;k--) {
	result[k]=residuo[k+2];
	for (j=k+1;j>=k;j--){ residuo[j] -= result[k]*cuadratic[j-k];}
	}
	for (j=2;j<=n;j++) residuo[j]=T(0);

	return result;


}//end deflate2

template<typename T>
boost::math::tools::polynomial<T>  deflate2(const boost::math::tools::polynomial<T> & poly, const
		std::complex<T>& root,boost::math::tools::polynomial<T> & residuo){



	T c=(root.real()*root.real())+(root.imag()*root.imag());

	T b= (T)((double)-2*root.real());
	T a=1;
	boost::math::tools::polynomial<T>  result=poly;
	boost::math::tools::polynomial<T>  cuadratic={{c,b,a}};
	int n=poly.size()-1;
	residuo=poly;

	int k,j;
	for (j=0;j<=n;j++) {

	result[j]=T(0);
	}

	for (k=n-2;k>=0;k--) {
	result[k]=residuo[k+2];
	for (j=k+1;j>=k;j--){ residuo[j] -= result[k]*cuadratic[j-k];}
	}
	for (j=2;j<=n;j++) residuo[j]=T(0);

	return result;


}//end deflate2



#endif /* DEFLATEALGORITHMS_H_ */
