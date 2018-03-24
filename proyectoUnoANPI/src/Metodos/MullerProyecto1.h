/*
 * MullerProyecto1.h
 *
 *  Created on: 22 mar. 2018
 *      Author: acacia
 */

#ifndef METODOS_MULLERPROYECTO1_H_
#define METODOS_MULLERPROYECTO1_H_


#include <iostream>
#include <complex>
#include "math.h"
#include <boost/math/tools/polynomial.hpp>
#include "../deflateAlgorithms.h"
using namespace boost::math::tools;
using namespace boost::math;

template <typename T>
T mullerMetod (const polynomial<T>& poly,T xr,bool polish,std::complex<T> cmplxRoot){
	if(poly.degree()<=1){
		return 0;
	}
	 float h = 1;
	    float eps = 0.001;
	    int maxit = 1000;
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
	    while(1){
	        iter++;
	        h0 = x1 - x0;
	        h1 = x2 - x1;
	        d0 = (poly.evaluate(x1) - poly.evaluate(x0))/h0;
	        d1 = (poly.evaluate(x2) - poly.evaluate(x1))/h1;
	        a = (d1 - d0)/ (h1 + h0);
	        b = a*h1 + d1;
	        c = poly.evaluate(x2);
	        rad = sqrt(b*b -4*a*c);
	        if (isnan(rad)){
	        		std::cout << "Tenemos una raíz compleja" << std::endl;
	        		std::complex<T> cRad(0,sqrt(fabs(b*b -4*a*c)));
	        		std::complex<T> cDen(0,0);
	        		std::complex<T> cDxr(0,0);
	        		if (std::abs(b+cRad) > std::abs(b-cRad)){
	        				cDen = b+cRad;
	        			}
	        		else{
	        				cDen = b-cRad;
	        			}
	        		cDxr = -2*c / cDen;
	        		cmplxRoot = x2 + cDxr;
	        		return NAN;
	        	}
	        if (fabs(b+rad) > fabs(b-rad)){
	            den = b+rad;
	        }
	        else {
	            den = b-rad;
	        }
	        dxr = -2*c / den;
	        xr = x2 + dxr;

	        if(polish){
	        	xr=mullerMetod<T>(poly,xr,false,cmplxRoot);
	        }

	        return xr;
	    }
	}


template<typename T>
std::vector<T> deflatedMuller(polynomial<T>& poly,T xr,bool polish){
	int polyGrade=poly.degree(); //input polynomial grade
	std::vector<T> roots(polyGrade); //output vector of the calculated roots
	T resi =0;
	polynomial<T> residuo={{0,0}};
	std::complex<T> complx(0,0);
	for(int i=0;i<polyGrade;i++){
		roots[i]=mullerMetod(poly, xr, polish,complx);
		if(isnan(roots[i])){
			poly=deflate2<T>(poly,complx,residuo);

			i++;
			roots[i]=NAN;
		}
		else{
		poly=deflate<T>(poly,roots[i],resi);

		}
	}
	std::cout<<"Raíces reales encontradas usando muller:"<<std::endl;
	for(int i=0;i<(int)roots.size();i++){

		if(!isnan(roots[i])){
		std::cout<<roots[i]<<std::endl;
		}
	}
	return roots ;
}


#endif /* METODOS_MULLERPROYECTO1_H_ */
