#include <iostream.h>
#include <math.h>
#include <unistd.h>
#include "lisasim.h"
#include "lisatest.h"

// this program compares the incoming Newtonian signal generated by the LISA Simulator
// with a comparable signal, generated by Synthetic LISA

int main(int argc, char **argv) {
    // find orbital separation
    
    double forb = 0.5*fgw;
    double R = pow( G*(M1+M2)*pow(1.0/(2.0*M_PI*forb),2) , 1.0/3.0);
    
    // find overall amplitude
    
    double A = 2.0*M1*M2/(r*R) * pow(G/(c*c),2);
    
    // I don't have a minus in my definition of inclination as in Newtonian.c, so I pass pi - inc [cos(pi-inc) = -cos(inc)]
    // need to adjust my frequency for the different definition of year in the Montana and JPL codes
    // need to convert my ecliptic latitude to the Montana latitude (which is Pi/2 - ...)
    // need to use minus my polarization angle
        
    SimpleBinary mybinary(fgw*1.000702365550482, phase, (M_PI-inc), A, 0.5*M_PI - theta, phi, -psi);
    
    ofstream myout("newtonian.txt");
    
    // remember that the LISA functions take times as fractions of a year
  
    long samples = (long)pow(2.0,nsource);
    
    // set timestep (in seconds) for 1 year of observations
    // and extend signal as per the LISA simulator
    
    double Dt = ObsTime / samples;
    long Ntot = samples + 2*(int)ceil(2.*AU/(c*Dt));
    double toff = -Dt*ceil(2.*AU/(c*Dt));

    // compute and write the signals at the Sun
 
    cout << "Computing " << Ntot << " samples" << endl;
 
    // if we use only hp and hc we need to do the polarization mixing explicitly
 
    double c2p = cos(2*psi);
    double s2p = sin(2*psi);
 
    for(int i=0;i<Ntot;i++) {
        double time = (toff + i*Dt)/year;
        
        double mhp = mybinary.hp(time);
        double mhc = mybinary.hc(time);
        
        myout << c2p*mhp + s2p*mhc << " " << -s2p*mhp + c2p*mhc << endl;
        
        if (i == (int)floor(Ntot/5.))    printf("  20 percent done\n");
        if (i == (int)floor(2.*Ntot/5.)) printf("  40 percent done\n");
        if (i == (int)floor(3.*Ntot/5.)) printf("  60 percent done\n");
        if (i == (int)floor(4.*Ntot/5.)) printf("  80 percent done\n");
        if (i == Ntot-1)                 printf(" 100 percent done\n");
    }
}
