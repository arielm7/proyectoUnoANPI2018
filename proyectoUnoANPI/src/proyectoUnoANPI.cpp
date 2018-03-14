//============================================================================
// Name        : proyectoUnoANPI.cpp
// Author      : Ariel
// Version     :
// Copyright   : Your copyright notice
// Description : Main
//============================================================================


#include <iostream>
#include <boost/math/tools/polynomial.hpp>
#include <iomanip>
#include <complex>
#include <cmath>
#include "deflateAlgorithms.cpp"
using namespace std;

using namespace boost::math;
using namespace boost::math::tools; // for polynomial
int main() {
	complex<double> x=7;
	complex<double> y=9;
	polynomial<double> const c{{4,0,1}};
	polynomial<double> resi ={{0,0,0}};
	boost::array<double, 3> a = {{-4, 0,1}};

	complex<double> const j(1,-3);
	cout << deflate2(c,j,resi) << endl; // prints hola
	cout<<"residuo: "<<resi<<endl;
	return 0;
}
