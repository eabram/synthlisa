//
//  August 2018
//  Ernst-Jan Buis
// 
//   $Id: 
//

#include <fstream>
#include <iostream>
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

//#include "lisasim-lisa.h"
//#include "lisasim-tens.h"
//#include "lisasim-wave.h"
//#include "lisasim-tdinoise.h"
//#include "lisasim-tdisignal.h"

//#include "SimpleBinary.h"


//#include <gsl/gsl_const_mksa.h>

int main(int argc, const char* argv[])
{

  const double c = 300000000;//GSL_CONST_MKSA_SPEED_OF_LIGHT;
  const double armlength = 2.5e9/c;

  // read orbit file
  std::vector<double> sc0, sc1, sc2;
  double sc0_x, sc0_y, sc0_z;
  double sc1_x, sc1_y, sc1_z;
  double sc2_x, sc2_y, sc2_z;
  double julian_date, delta_t, old_date = 0;

  std::ifstream infile("positions.txt");

  while(1)
    {
      infile >> julian_date >> sc0_x >> sc0_y >> sc0_z >> sc1_x >> sc1_y >> sc1_z >> sc2_x >> sc2_y >> sc2_z;
      if (!infile.good()) break;
      sc0.push_back(sc0_x/c);
      sc0.push_back(sc0_y/c);
      sc0.push_back(sc0_z/c);

      sc1.push_back(sc1_x/c);
      sc1.push_back(sc1_y/c);
      sc1.push_back(sc1_z/c);

      sc2.push_back(sc2_x/c);
      sc2.push_back(sc2_y/c);
      sc2.push_back(sc2_z/c);
      delta_t = julian_date - old_date;
      old_date = julian_date;
      //      std::cout << julian_date << "\t" << sc0_x << "\t" << delta_t << std::endl;

    }

  delta_t *= 24 * 3600.;
  infile.close();
  long length = sc1.size();

  double *pos0 = &sc0[0]; // convert std vector to C++ array
  double *pos1 = &sc1[0];
  double *pos2 = &sc2[0];
  
  SampledLISA *lisamission = new SampledLISA(pos0, length, pos1, length, pos2, length, delta_t, 0, 1);
  double st_laser_noise[6] = {1,1,1,1,1,1};
  double st_shot_noise[6]  = {1,1,1,1,1,1};
  double st_pm_noise[6]    = {1,1,1,1,1,1};

  double sd_laser_noise[6] = {0,0,0,0,0,0};
  double sd_shot_noise[6]  = {1e-46, 1e-46, 1e-46, 1e-46, 1e-46, 1e-46};
  double sd_pm_noise[6]    = {0,0,0,0,0,0};
  
  /*TDIpointing *pointingTDI = new TDIpointing(lisamission,
  					     st_pm_noise,    sd_pm_noise,
  					     st_shot_noise,  sd_shot_noise,
  					     st_laser_noise, sd_laser_noise);

  
  SimpleBinary *mysystem = new SimpleBinary( 1.0e-4,
					     0.0,
					     0,
					     1.4,
					     0,
					     0,
					     0);

    // create TDInoise object
  TDIsignal *eccentricTDI = new TDIsignal(lisamission, mysystem);

  
  for (int i =0; i <  3e8; i+= 1e4)
    {
      Vector p1,p2;
      Vector u,v;
      lisamission->putp(p1, 1, i);
      lisamission->putp(p2, 2, i);

      lisamission->putn(v, 1, i);
      lisamission->putn(u, 2, i);
      double angle = dotproduct(u, v);

      Vector uu,vv;
      lisamission->putn(vv, 1, i);
      lisamission->putn(uu, 3, i);
      double angle2 = dotproduct(uu, vv);

      Vector uuu,vvv;
      lisamission->putn(vvv, 2, i);
      lisamission->putn(uuu, 3, i);
      double angle3 = dotproduct(uuu, vvv);

      std::cout
 	<< i
	<< "\t" << angle
	<< "\t" << angle2
	<< "\t" << angle3
	<< "\t" << lisamission->armlength(1, i)
	<< "\t" << p1[0]
	<< "\t" << p1[1]
	<< "\t" << p1[2]
//	<< "\t" << p2[0]
//	<< "\t" << p2[1]
//	<< "\t" << p2[2]
	// << "\t" << eccentricTDI->Xm(i)
	// << "\t" << eccentricTDI->Ym(i)
	// << "\t" << eccentricTDI->Zm(i)
// 	      << "\t" << pointingTDI->y123(i)
// 	      << "\t" << pointingTDI->y231(i)
// 	      << "\t" << pointingTDI->y312(i)
// 	      << "\t" << pointingTDI->z123(i)
// 	      << "\t" << pointingTDI->z213(i)
//	      << "\t" << pointingTDI->z312(i)
	      << std::endl;
    }
      

 */  
}
