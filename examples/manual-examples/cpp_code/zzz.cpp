//Test programm
#include <iostream>
using namespace std;


#include <fstream>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "/home/ester/git/synthlisa/lisasim/lisasim-lisa.h"
#include "/home/ester/git/synthlisa/lisasim/lisasim-tens.h"
#include "/home/ester/git/synthlisa/lisasim/lisasim-wave.h"
#include "/home/ester/git/synthlisa/lisasim/lisasim-tdinoise.h"
#include "/home/ester/git/synthlisa/lisasim/lisasim-tdisignal.h"

//#include "/home/ester/git/synthlisa/lisasim/SimpleBinary.h"


//#include <gsl/gsl_const_mksa.h>



int main()
{
	//const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
  	const double c = 300000000;
	const double armlength = 2.5e9/c;
	const double scale = 1000; //...adjusted;
	const double year2sec = 32536000;

	std::vector<double> sc1, sc2, sc3,t_vec;
	double p1x, p1y, p1z;
	double p2x, p2y, p2z;
	double p3x, p3y, p3z;
	double t, dt;

	std::ifstream infile("positions.txt");
	int count=0;

	while(1){
		//std::cout << "Reading";	
		infile >> t >> p1x >> p1y >> p1z >> p2x >> p2y >> p2z >> p3x >> p3y >> p3z;
		if (!infile.good()) break;
		count = count +1;
		sc1.push_back(p1x/c);
		sc1.push_back(p1y/c);
		sc1.push_back(p1z/c);
		
		sc2.push_back(p2x/c);
		sc2.push_back(p2y/c);
		sc2.push_back(p2z/c);

		sc3.push_back(p3x/c);
		sc3.push_back(p3y/c);
		sc3.push_back(p3z/c);
	
		t_vec.push_back(t*24*3600);
		
		//std::cout << p1x << std::endl;
		
	}
	std::cout << sc1[0] << endl;
	std::cout << "DONE" << endl;
	infile.close();
	dt = t_vec[1]-t_vec[0]; //24*3600;
	std::cout<<"dt = "<<dt<<endl;
  	long length = sc1.size()/3;

	double *p1 = &sc1[0]; // convert std vector to C++ array
        double *p2 = &sc2[0];
        double *p3 = &sc3[0];
	
	double pos1[length][3];	
	double pos2[length][3];
	double pos3[length][3];


	

	for (int i = 0; i<length ; i++){
		//std::cout<<sc1[i]<< sc1[i+1] <<  sc1[i+2]<<std::endl;
		//std::cout << i <<endl;
		for (int j=0; j<3; j++){
			//std::cout << i;
			pos1[i][j]=*(p1+i*3+j);
			pos2[i][j]=*(p2+i*3+j);
			pos3[i][j]=*(p3+i*3+j);


			//	sc1[i*3+j];
			//pos2[i][j]=sc2[i*3+j];
			//pos3[i][j]=sc3[i*3+j];
		}
		//std::cout<<&pos1[i]<<endl;
	}
	
	std::cout<< *p1 << endl;
	long l = sc1.size();

  	//SampledLISA *lisamission = new SampledLISA(p1,l,p2,l,p3,l, dt, 0, 1);	
  	SampledLISA *lisamission = new SampledLISA(*pos1,length,*pos2,length,*pos3,length, dt, 0, 1);
	double L = 16;
	//OriginalLISA *lisaoriginal = new OriginalLISA(L,L,L);
	
	
	
	
	return 0;
}




