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
//#include "deflateAlgorithms.cpp"
#include "MullerProyecto1.cpp"
using namespace std;

using namespace boost::math;
using namespace boost::math::tools; // for polynomial
int main() {

	polynomial<double> const c{{-6,1,1}};
	double resi =0;
	double const j=2.01;
	//cout << deflate(c,j,resi) << endl; // prints hola
	//cout<<"residuo: "<<resi<<endl;

	double point = -5.3;
	polynomial<double> poly{{-10,-11,4,1}};
	std::vector<double> raices = deflatedMuller<double>(poly, point);
	std::cout <<"raices: "<< raices[1] << std::endl;
	return 0;
}
