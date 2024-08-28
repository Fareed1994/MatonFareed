#include "RGEsRunner.h"
#include "RGEsDefiner.h"
#include "RungeKutta.h"
#include <Minuit2/MnMigrad.h>
#include "ChiSquared.h"

namespace RGEsGlobals {
double penaltyFactor = 1.0;
}

VectorXcd RunRGEs(const std::vector<double>& parameters) {
    std::cout<<"Running!"<<std::endl;
    double pi = M_PI;
    double gGUT, MGUT, MEW, Lambda, TanBeta, M12, A, A0, M0;
    dcomplex MHu2, MHd2, mu;
    Matrix3cd Yu, Yd, Ye, TYu, TYd, TYe, Au, Ad, Ae, Mq2, Ml2, Mu2, Md2, Me2;
    Vector3cd M;
    Vector3d g;
    double PlankScale = 1e18; //GeV
    
    double b1 = 31./5;
    double b2 = 1.;
    double b3 = -3.;
    double M1, M2, M3;
    
    double e3, lam, re, ret, rep, rxi, rdel, phidel, ralp, phialp, rbet,
           phibet, phie, phiet, phiep, phixi, rth, phith;
           

    std::complex<double> del, alp, bet, e, et, ep, xi, th;

    // Mapping the parameters into the model:
    gGUT = sqrt(4.0 * pi / parameters[0]);
    MGUT = parameters[1] * 1.e+16;
    MEW = parameters[2];
    Lambda  = parameters[3];
    TanBeta = parameters[4];
    mu      = parameters[5];
    M12     = parameters[6];
    A       = parameters[7];
    M0      = parameters[8];
    MHd2    = dcomplex(parameters[9] * M0 * M0, 0.);
    MHu2    = dcomplex(parameters[10] * M0 * M0, 0.);
    A0      = parameters[11];
    re      = parameters[13];
    ret     = parameters[14];
    rep     = parameters[15];
    rxi     = parameters[16];
    rdel    = parameters[17];
    phidel  = parameters[18];
    ralp    = parameters[19];
    phialp  = parameters[20];
    rbet    = parameters[21];
    phibet  = parameters[22];
    phie    = parameters[23];
    phiet   = parameters[24];
    phiep   = parameters[25];
    phixi   = parameters[26];
    rth     = parameters[27];
    phith   = parameters[28];
    
    g << gGUT, gGUT, std::sqrt(1 + parameters[12]) * gGUT;

    M1 = (1+(gGUT*gGUT*b1*A)/(16*pi*pi)*std::log(PlankScale/M0))*M12;
    M2 = (1+(gGUT*gGUT*b2*A)/(16*pi*pi)*std::log(PlankScale/M0))*M12;
    M3 = (1+(gGUT*gGUT*b3*A)/(16*pi*pi)*std::log(PlankScale/M0))*M12;

    M << dcomplex(M1, 0), dcomplex(M2, 0), dcomplex(M3, 0);
    
    
    del = std::complex<double> (rdel*cos(phidel), rdel*sin(phidel));
    alp = std::complex<double> (ralp*cos(phialp), ralp*sin(phialp));
    bet = std::complex<double> (rbet*cos(phibet), rbet*sin(phibet));
    e   = std::complex<double> (re*cos(phie), re*sin(phie));
    et  = std::complex<double> (ret*cos(phiet), ret*sin(phiet));
    ep  = std::complex<double> (rep*cos(phiep), rep*sin(phiep));
    xi  = std::complex<double> (rxi*cos(phixi), rxi*sin(phixi));
    th  = std::complex<double> (rth*cos(phith), rth*sin(phith));
    
    Yu << 0, ep * ( 1./(1.+alp-(4./3.)*bet)-1./(1.+alp+(1./3.)*bet) )-th*(          1./(1.+alp-(4./3.)*bet)+1./(1.+alp+(1./3.)*bet) ), -xi / (1.+alp-(4./3.)*bet),
        -ep * ( 1./(1.+alp-(4./3.)*bet)-1./(1.+alp+(1./3.)*bet) )-th*( 1./(1.+alp-(4./3.)*bet)+1./(1.+alp+(1./3.)*bet) ), et * ( 1./(1.+alp-(4./3.)*bet)-1./(1.+alp+(1./3.)*bet) ), -e / (1.+alp-(4./3.)*bet),
        xi /(1.+alp+(1./3.)*bet), e /(1.+alp+(1./3.)*bet), dcomplex(Lambda, 0);

    Yd << 0, ep * ( 1./(1.-3.*alp+(2./3.)*bet)-1./(1.+alp+(1./3.)*bet) )-th*( 1./(1.-3.*alp+(2./3.)*bet)+1./(1.+alp+(1./3.)*bet) ), -xi / (1.-3.*alp+(2./3.)*bet),
        -ep * ( 1./(1.-3.*alp+(2./3.)*bet)-1./(1.+alp+(1./3.)*bet) )-th*( 1./(1.-3.*alp+(2./3.)*bet)+1./(1.+alp+(1./3.)*bet) ), et * ( 1./(1.-3.*alp+(2./3.)*bet)-1./(1.+alp+(1./3.)*bet) ), -e / (1.-3.*alp+(2./3.)*bet),
        xi / (1.+alp+(1./3.)*bet), e /(1.+alp+(1./3.)*bet), dcomplex(Lambda, 0);

    Ye << 0, -ep * ( 1./(1.-3.*alp-bet)-1./(1.+alp+2.*bet) )-th*( 1./(1.-3.*alp-bet)+1./(1.+alp+2.*bet) ), 3. * xi  / (1.+alp+2.*bet),
        ep * ( 1./(1.-3.*alp-bet)-1./(1.+alp+2.*bet) )-th*th*( 1./(1.-3.*alp-bet)+1./(1.+alp+2.*bet) ), 3. * et * ( 1./(1.-3.*alp-bet)-1./(1.+alp+2.*bet) ), 3. * e / (1.+alp+2.*bet),
        -3.* xi  / (1.-3.*alp-bet), -3.* e / (1.-3.*alp-1.*bet), dcomplex(Lambda, 0);

    TYu = A0 * Yu;
    TYd = A0 * Yd;
    TYe = A0 * Ye;
    
    Yu = Yu.transpose().eval();     // In model.cpp, Yukawas are defined with doublets on the left,
    Yd = Yd.transpose().eval();     // but RGEs are written with doublets on the right! We Transpose
    Ye = Ye.transpose().eval();     // them here. Furthermore, Yukawas entering RGEs should be complex
                            // conjugate but I think it doesn't make any difference, as far as
    TYu = TYu.transpose().eval();   // we define diagonalization  matrices and VCKM properly, and as far
    TYd = TYd.transpose().eval();   // as soft susy breaking pars are real.
    TYe = TYe.transpose().eval();

    Mq2 << dcomplex(std::pow(M0, 2), 0), 0, 0,
           0, dcomplex(std::pow(M0, 2), 0), 0,
           0, 0, dcomplex(std::pow(M0, 2), 0);

    Ml2 << dcomplex(std::pow(M0, 2), 0), 0, 0,
           0, dcomplex(std::pow(M0, 2), 0), 0,
           0, 0, dcomplex(std::pow(M0, 2), 0);

    Mu2 << dcomplex(std::pow(M0, 2), 0), 0, 0,
           0, dcomplex(std::pow(M0, 2), 0), 0,
           0, 0, dcomplex(std::pow(M0, 2), 0);

    Md2 << dcomplex(std::pow(M0, 2), 0), 0, 0,
           0, dcomplex(std::pow(M0, 2), 0), 0,
           0, 0, dcomplex(std::pow(M0, 2), 0);

    Me2 << dcomplex(std::pow(M0, 2), 0), 0, 0,
           0, dcomplex(std::pow(M0, 2), 0), 0,
           0, 0, dcomplex(std::pow(M0, 2), 0);

    // Mapping the parameters into a vector y(t):
    VectorXcd ystart(107);

    int k = 0;
    for(int i=0; i<=2; i++)
    {
        for (int j=0; j<=2; j++)
        {
            ystart(k) = Yu(i,j);
            ystart(9+k) = Yd(i,j);
            ystart(18+k) = Ye(i,j);
            ystart(27+k) = TYu(i,j);
            ystart(36+k) = TYd(i,j);
            ystart(45+k) = TYe(i,j);
            ystart(54+k) = Mq2(i,j);
            ystart(63+k) = Ml2(i,j);
            ystart(72+k) = Mu2(i,j);
            ystart(81+k) = Md2(i,j);
            ystart(90+k) = Me2(i, j);
            k++;
        }
        ystart(99+i) = M(i);
        ystart(102+i) = std::real(g(i));
    }
    ystart(105) = MHu2;
    ystart(106) = MHd2;
    //ystart(107) = Bmu;

    //parameters for integration:
    double eps = 1.0e-6;           // Error tolerance
    double hmin = 0.0;             // Minimum step size
    double h1 = (std::log(MGUT) - std::log(MEW))/1.0e4; // Initial step size
    // Integrate
    int status = integrateOdes(ystart, std::log(MGUT), std::log(MEW), eps, h1, hmin, RGEsystem);
    if (std::abs(status-1.0) > eps) {
        RGEsGlobals::penaltyFactor = 1.0;
        std::cout << "--- Integration iteration successful! ---" << std::endl;
        return ystart;
    } else {
        RGEsGlobals::penaltyFactor = 1.e6;
        std::cout << "xxx   Integration iteration failed.   xxx" << std::endl;
        return ystart;
    }
}







	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	 
	   
