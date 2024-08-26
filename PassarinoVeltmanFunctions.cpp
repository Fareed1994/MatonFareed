//--------------------------------------
// Passarino - Veltman Functions
// Pierce et al, T. Blazek - notations
//--------------------------------------
#include "PassarinoVeltmanFunctions.h"
// ****************************************************************************
// A0 function 
// T. Blazek notation, opposite sign in Pierce et al

double A0(double M)
{
  double Q;
  Q = 91.188;
//#K
//  Q = 100.00;

  if (M == 0.) return 0.;
  else  return - M*M*(1. - 2.*log(M/Q));
}

// ****************************************************************************
// B0 function
// T. Blazek notation =  Pierce et al

double B0(double p, double m1, double m2)
{

//  std::cout<<"m1 = "<<m1<<"    m2 = "<<m2<<"     p = "<<p<<std::endl;

  double Q;
  Q = 91.188;
//#K
//Q = 100.0;

  if (m1 == 0.)
  {
    if (p == 0. && m2 == 0.) return 1.; // B0 is infinite, but it doesn't
                                        // matter since B0 is multiplied by 0!
    else if (p == 0.)        return 1. - 2.*log(m2/Q);
    else if (p == m2)        return 2. - 2.*log(m2/Q);
    else if (m2 == 0.)       return 2. - 2.*log(p/Q);
    else                     return 2. - 2.*log(m2/Q) +
                             (1. - m2*m2/(p*p)) * log(m2*m2/fabs(m2*m2-p*p));
  }
  
  else if (m2 == 0.)
  {
    if      (p == 0.)  return 1. - 2.*log(m1/Q);
    else if (p == m1)  return 2. - 2.*log(m1/Q);
    else return 2. - 2.*log(m1/Q) +
                (1. - m1*m1/(p*p)) * log(m1*m1/fabs(m1*m1-p*p));
  }

  else if (m1 == m2)
  {
    if (p == 0.)        return -2.*log(m1/Q);
    else if (p > 2.*m1) return -2.*log(m1/Q) + 2. - (1./p)*sqrt(p*p - 4.*m1*m1)
                        * log((p+sqrt(p*p - 4*m1*m1))/(p-sqrt(p*p - 4*m1*m1)));
    else                return -2.*log(m1/Q) + 2. - (2./p)*sqrt(4.*m1*m1 - p*p)
                        * atan(p/sqrt(4.*m1*m1 - p*p));
  }

  else if (p == 0.)
    return - log(m1/Q) - log(m2/Q) + 1.
           - (m1*m1 + m2*m2) * log(m1/m2) / ((m1-m2)*(m1+m2));

  else if (p*p < (m1-m2)*(m1-m2))
    return - log(m1/Q) - log(m2/Q) + 1.
           - (m1*m1 + m2*m2) * log(m1/m2) / ((m1-m2)*(m1+m2))
           + 1. + ((m1-m2)*(m1+m2)/(p*p) - (m1*m1 + m2*m2)/((m1-m2)*(m1+m2)))
           * log(m2/m1)
           + sqrt(((m1+m2)*(m1+m2) - p*p)*((m1-m2)*(m1-m2) - p*p))
           * log( (sqrt((m1+m2)*(m1+m2) - p*p) + sqrt((m1-m2)*(m1-m2) - p*p))
           / (sqrt((m1+m2)*(m1+m2) - p*p) - sqrt((m1-m2)*(m1-m2) - p*p)) )
           / (p*p);

  else if ( p*p > (m1-m2)*(m1-m2) && p*p < (m1+m2)*(m1+m2) )
    return - log(m1/Q) - log(m2/Q) + 1.
           - (m1*m1 + m2*m2) * log(m1/m2) / ((m1-m2)*(m1+m2))
           + 1. + ((m1-m2)*(m1+m2)/(p*p) - (m1*m1 + m2*m2)/((m1-m2)*(m1+m2)))
           * log(m2/m1)
           - sqrt(((m1+m2)*(m1+m2) - p*p)*(p*p - (m1-m2)*(m1-m2)))
           * atan( sqrt(p*p - (m1-m2)*(m1-m2))/sqrt((m1+m2)*(m1+m2) - p*p))
           * 2. / (p*p);

  else if ( p*p > (m1+m2)*(m1+m2) )
    {
//std::cout<<"sqrt(p*p - (m1-m2)*(m1-m2)) = "<<sqrt(p*p - (m1-m2)*(m1-m2))<<"\n";
//std::cout<<"sqrt(p*p - (m1+m2)*(m1+m2)) = "<<sqrt(p*p - (m1+m2)*(m1+m2))<<"\n";
    return - log(m1/Q) - log(m2/Q) + 1.
           - (m1*m1 + m2*m2) * log(m1/m2) / ((m1-m2)*(m1+m2))
           + 1. + ((m1-m2)*(m1+m2)/(p*p) - (m1*m1 + m2*m2)/((m1-m2)*(m1+m2)))
           * log(m2/m1)
           - sqrt((p*p - (m1+m2)*(m1+m2))*(p*p - (m1-m2)*(m1-m2)))
           * log( (sqrt(p*p - (m1+m2)*(m1+m2)) + sqrt(p*p - (m1-m2)*(m1-m2)))
           / (sqrt(p*p - (m1-m2)*(m1-m2)) - sqrt(p*p - (m1+m2)*(m1+m2))) )
           / (p*p);

    }

  else {
         //std::cout<<"No case in B0(p, m1, m2)!"<<std::endl;
         RGEsGlobals::penaltyFactor *= 1.e5;
         return 0;
       }
}

// ****************************************************************************8
// B1 function
// T. Blazek notation, opposite sign in Pierce et al
// doesn't work if p = 0 and one of the masses is 0

    
double B1(double p, double M, double m)
{
  double x, Q;
  Q = 91.188;
//#K
//Q = 100.0;

  if (M == m) return -0.5*B0(p,M,m);

  else if (p != 0.) return (A0(m) - A0(M) + (m*m-M*M-p*p)*B0(p,M,m))/(2.*p*p);

  else if (M == 0. && m == 0.) return 1.; // B0 is infinite, but it doesn't
                                        // matter since B0 is multiplied by 0!
  else if (M == 0.) return -0.25 + log(m/Q);

  else if (m == 0.) return -0.75 + log(M/Q);

  else if (M > m)
  {
    x = m*m/(M*M);
    return log(M/Q) - .25 - .5/(1.-x) - .5*log(x)/((1.-x)*(1.-x)) + .5*log(x);
  }
  else
  {
    x = m*m/(M*M);
    return log(m/Q) - .25 - .5/(1.-x) - .5*log(x)/((1.-x)*(1.-x));
  }
 }


// ****************************************************************************8
// B22 and tildeB22 functions
// T. Blazek notation, opposite sign in Pierce et al

double B22(double p, double m1, double m2)
{
  if (p == 0.) return ( 0.5*A0(m1) + 0.5*A0(m2) - m1*m1 - m2*m2 -
               (1.5*m1*m1 + 0.5*m2*m2) * B0(p,m1,m2)
               + (m2*m2 - m1*m1) * B1(p,m1,m2) )/6.;

  else     return ( 0.5*A0(m1) + 0.5*A0(m2) + p*p/3. - m1*m1 - m2*m2 +
           (0.5*p*p - m1*m1 - m2*m2 + (m1*m1 - m2*m2)*(m1*m1 - m2*m2)/(2.*p*p))
           * B0(p,m1,m2)
           +(m1*m1 - m2*m2)*(A0(m1) - A0(m2))/(2.*p*p) )/6.;
}

double tB22(double p, double m1, double m2)
{
//  std::cout<<"m1 = "<<m1<<"    m2 = "<<m2<<"     p = "<<p<<std::endl;
  return B22(p, m1, m2) -  0.25*A0(m1) - 0.25*A0(m2);
}


// ****************************************************************************8
// Other functions
// Pierce et al with appropriate signs!!!!!

double F(double p, double m1, double m2)
{
  return - A0(m1) + 2.*A0(m2) - (2.*p*p + 2.*m1*m1 - m2*m2) * B0(p,m1,m2);
}

double G(double p, double m1, double m2)
{
  return (p*p - m1*m1 - m2*m2) * B0(p,m1,m2) + A0(m1) + A0(m2);
}

double H(double p, double m1, double m2)
{
  return -4.*B22(p,m1,m2) + G(p,m1,m2);
}


double C0(double m1, double m2, double m3)
{
  return ( 2.*m2*m2*log(m2/m1)/((m1-m2)*(m1+m2)) -
           2.*m3*m3*log(m3/m1)/((m1-m3)*(m1+m3)) ) / ((m2-m3)*(m2+m3));
}

double D0(double m1, double m2, double m3, double m4)
{
  return ( C0(m1,m3,m4) - C0(m2,m3,m4)  ) / ((m1-m2)*(m1+m2));
}
