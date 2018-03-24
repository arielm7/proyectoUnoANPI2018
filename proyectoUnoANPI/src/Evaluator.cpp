
#include <iostream>
#include <complex>
#include "math.h"
#include <boost/math/tools/polynomial.hpp>
#include "Metodos/MullerProyecto1.h"
using namespace boost::math::tools;
using namespace boost::math;

void evaluateMethod(){
	std::cout << "Se evaluará el siguiente polinomio: 6*x^3+561.043x^2+1884.15x-871.454" << std::endl;	
	std::cout << "****Sin pulir las raices****" << std::endl;	
	int valorR1 = 1;
	int valorR2 = -2;
	double pointDouble = 3.0;
	polynomial<double> polyDouble{{-871.454,1884.15,561.043,6}};
	std::vector<double> vectorResults = deflatedMuller<double>(polyDouble, pointDouble, false);
	std::cout <<"Raíz 1 obtenida: "<< vectorResults.at(0)<<  "     " << "Raiz 1 real: 0.4118 " << std::endl; 
	std::cout <<"Raíz 2 obtenida: "<< vectorResults.at(1) << "     " << "Raiz 2 real: -3.9189" << std::endl; 
	std::cout <<"Raíz 3 obtenida: "<< vectorResults.at(2) << "     " << "Raiz 3 real: -90" <<std::endl; 
	std::cout <<"Errores" << std::endl; 
	std::cout << "Error raíz 1 = " << ((0.4118-vectorResults.at(0))/0.4118)*100 << std::endl;
	std::cout << "Error raíz 2 = " << ((-3.9189-vectorResults.at(1))/-3.9189)*100 << std::endl;
	std::cout << "Error raíz 3 = " << ((-90-vectorResults.at(2))/-90)*100 << std::endl;
	std::cout << "****Puliendo las raices****" << std::endl;
	polynomial<double> polyDoubleP{{-871.454,1884.15,561.043,6}};
	std::vector<double> vectorResultsPuliendo = deflatedMuller<double>(polyDoubleP, pointDouble, true);
	std::cout <<"Raíz 1 obtenida: "<< vectorResultsPuliendo.at(0)<<  "     " << "Raiz 1 real: 0.4118 " << std::endl; 
	std::cout <<"Raíz 2 obtenida: "<< vectorResultsPuliendo.at(1) << "     " << "Raiz 2 real: -3.9189" << std::endl; 
	std::cout <<"Raíz 3 obtenida: "<< vectorResultsPuliendo.at(2) << "     " << "Raiz 3 real: -90" <<std::endl; 
	std::cout <<"Errores" << std::endl; 
	std::cout << "Error raíz 1 = " << ((0.4118-vectorResultsPuliendo.at(0))/0.4118)*100 << std::endl;
	std::cout << "Error raíz 2 = " << ((-3.9189-vectorResultsPuliendo.at(1))/-3.9189)*100 << std::endl;
	std::cout << "Error raíz 3 = " << ((-90-vectorResultsPuliendo.at(2))/-90)*100 << std::endl;
	std::cout << "****Utilizando floats****" << std::endl;
	polynomial<float> polyDoubleF{{-871.454,1884.15,561.043,6}};
	std::vector<float> vectorResultsFloats = deflatedMuller<float>(polyDoubleF, pointDouble, true);
	std::cout <<"Raíz 1 obtenida: "<< vectorResultsFloats.at(0)<<  "     " << "Raiz 1 real: 0.4118 " << std::endl; 
	std::cout <<"Raíz 2 obtenida: "<< vectorResultsFloats.at(1) << "     " << "Raiz 2 real: -3.9189" << std::endl; 
	std::cout <<"Raíz 3 obtenida: "<< vectorResultsFloats.at(2) << "     " << "Raiz 3 real: -90" <<std::endl; 
	std::cout <<"Errores" << std::endl; 
	std::cout << "Error raíz 1 = " << ((0.4118-vectorResultsFloats.at(0))/0.4118)*100 << std::endl;
	std::cout << "Error raíz 2 = " << ((-3.9189-vectorResultsFloats.at(1))/-3.9189)*100 << std::endl;
	std::cout << "Error raíz 3 = " << ((-90-vectorResultsFloats.at(2))/-90)*100 << std::endl;
	
	
	
	//std::cout << "Usando precisión float" << std::endl;	
	//float pointFloat = 1.0;
	//polynomial<float> const polyFloat{{-18,9,7,1,1}};
	//mullerMetod(polyFloat, pointFloat, 1); 
	//std::cout << "Sin pulir las ráices" << std::endl;	
	//float pointSinPulir = 1.0;
	//polynomial<float> const polySinPulir{{-18,9,7,1,1}};
	//mullerMetod(polySinPulir, pointSinPulir, 0); 
	//std::cout << "Puliendo raíces" << std::endl;	
	//float pointPuliendo = 1.0;
	//polynomial<float> const polyPuliendo{{-18,9,7,1,1}};
	//mullerMetod(polyPuliendo, pointPuliendo, 1); 
	
	
	
}
