/*
 * main.cpp
 *
 *  Created on: 12 mar. 2018
 *      Author: acacia
 */

#include <boost/math/tools/polynomial.hpp>

#include "Metodos/laguerre/deflatedLaguerre.h"
#include "Metodos/MullerProyecto1.h"
#include <stdlib.h>


int main(int argc, char *argv[]){
	int num=argc;
	cout << *argv[1] << endl;
	if(num<4){
		std::cout<<"faltan argumentos"<<std::endl;
		return 0;
	}
	else if(num==4){
		std::cout<<"se ingreso una constante no un polinomio"<<std::endl;
		return 0;
	}

	if(*argv[1]=='d'){
		double coef[argc-3];
		for(int i=3;i<argc;i++){
			cout << coef[i-3] << endl;
			coef[i-3]=std::atof(argv[i]);

		}
		boost::math::tools::polynomial<double>*  const poli= new boost::math::tools::polynomial<double>(coef, argc-4);
		deflatedLaguerre<double>* laguer = new deflatedLaguerre<double>();
		std::cout<<*poli<<std::endl;


		double x = 0;
		double point = 0;
		if(*argv[2]=='p'){
			laguer->initLaguerre(poli,x,true);
			deflatedMuller<double>(*poli, point,true);
		}
		else{
			laguer->initLaguerre(poli,x,false);
			deflatedMuller<double>(*poli, point,false);
		}





	}
	if(*argv[1]=='f'){
			double coef[argc-3];
			for(int i=3;i<argc;i++){
				cout << coef[i-3] << endl;
				coef[i-3]=std::atof(argv[i]);

			}
			boost::math::tools::polynomial<float>*  const poli= new boost::math::tools::polynomial<float>(coef, argc-4);
			deflatedLaguerre<float>* laguer = new deflatedLaguerre<float>();
			std::cout<<*poli<<std::endl;


			float x = 0;
			float point = 0;
			if(*argv[2]=='p'){
				laguer->initLaguerre(poli,x,true);
				deflatedMuller<float>(*poli, point,true);
			}
			else{
				laguer->initLaguerre(poli,x,false);
				deflatedMuller<float>(*poli, point,false);
			}





		}

	else if(*argv[1]=='c' && *argv[2]=='f'){
				double coef[argc-3];
				for(int i=3;i<argc;i++){
					cout << coef[i-3] << endl;
					coef[i-3]=std::atof(argv[i]);

				}
				boost::math::tools::polynomial<std::complex<float>>*  const poli= new boost::math::tools::polynomial<std::complex<float>>(coef, argc-4);
				deflatedLaguerre<std::complex<float>>* laguer = new deflatedLaguerre<std::complex<float>>();
				std::cout<<*poli<<std::endl;


				std::complex<float> x = 0;
				if(*argv[2]=='p'){
					laguer->initLaguerre(poli,x,true);
				}
				else{
					laguer->initLaguerre(poli,x,false);
				}





			}


	else if(*argv[1]=='c' && *argv[2]=='d'){
					double coef[argc-3];
					for(int i=3;i<argc;i++){
						coef[i-3]=std::atof(argv[i]);

					}
					boost::math::tools::polynomial<std::complex<double>>*  const poli= new boost::math::tools::polynomial<std::complex<double>>(coef, argc-4);
					deflatedLaguerre<std::complex<double>>* laguer = new deflatedLaguerre<std::complex<double>>();
					std::cout<<*poli<<std::endl;


					std::complex<double> x = 0;
					if(*argv[2]=='p'){
						laguer->initLaguerre(poli,x,true);
					}
					else{
						laguer->initLaguerre(poli,x,false);
					}





				}


	return 0;
}
