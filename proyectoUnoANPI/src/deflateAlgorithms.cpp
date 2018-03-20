/*
 * deflateAlgorithms.cpp
 *
 *  Created on: mar. 2018
 *      Author: ariel
 */

#include <boost/math/tools/polynomial.hpp>
#include <complex>
#include <iostream>

using namespace boost::math;
using namespace boost::math::tools; // for polynomial

template<typename T>
polynomial<T> deflate(const polynomial<T>& poly, const T& root, T& residuo){

	polynomial<T>  result=poly;
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
polynomial<T> deflate2(const polynomial<T>& poly, const
		std::complex<T>& root,polynomial<T>& residuo){

	polynomial<T> result=poly;
	polynomial<T> cuadratic={{(root.real()*root.real())+(root.imag()*root.imag()),-2*root.real(),1}};
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

