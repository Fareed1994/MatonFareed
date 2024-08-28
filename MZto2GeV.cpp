#include "MZto2GeV.hpp"

double beta0(int nf)
{
  return (11. - 2.*nf/3.)/4.;
}
  
double beta1(int nf)
{
  return (102. - 38.*nf/3.)/16.;
}
  
double beta2(int nf)
{
  return (2857./2. - 5033.*nf/18. + 325.*nf*nf/54.)/64.;
}

double bb1(int nf)
{
  return beta1(nf)/beta0(nf);
}

double bb2(int nf)
{
  return beta2(nf)/beta0(nf);
}

double gam0(int nf)
{
  return 1.;
}
 
double gam1(int nf)
{
  return (202./3. - 20.*nf/9.)/16.;
}
 
double gam2(int nf)
{
  return (1249. + (-2216./27. -  160.*ZETA_3/3.)*nf - 140.*nf*nf/81.)/64.;
}


double del(double x)
{
  const double pi = 3.141592653589793;
  return pi*pi*x/8. - 0.597*x*x+0.230*x*x*x;
}

void run_qcd(int nf, double a_3_i, double& a_3_f,
             double sc_i, double sc_f, double& R)
{

  const double pi = 3.141592653589793;
  double a, L, bet0, b1, b2, c0, c1, c2;
  double x, c_f, c_i;

  bet0 = beta0(nf);
  b1 = bb1(nf);
  b2 = bb2(nf);
  c0 = gam0(nf)/bet0;
  c1 = gam1(nf)/bet0;
  c2 = gam2(nf)/bet0;

  a = a_3_i/pi;
  L = 2.*std::log(sc_f/sc_i) +
      (1./a + b1*std::log(a) + (b2-b1*b1)*a - b1*b2*a*a)/bet0 +
      (b1/bet0) * std::log(bet0);

  if (L<=0.)
  {
    R = -1.;
    return;
  }

  a = 1./(bet0*L) - b1*std::log(L)/(bet0*bet0*L*L) +
      (b1*b1*(std::log(L)*std::log(L) - std::log(L) - 1) + b2)/(bet0*bet0*bet0*L*L*L) -
      3.*b1*b2*std::log(L)/(bet0*bet0*bet0*bet0*L*L*L*L);

  a_3_f = a * pi;

  x = a_3_f/pi;
  c_f = (1. + (c1-b1*c0)*x +
         0.5*((c1-b1*c0)*(c1-b1*c0) + c2 - b1*c1 + b1*b1*c0 - b2*c0)*x*x +
         ((c1-b1*c0)*(c1-b1*c0)*(c1-b1*c0)/6. +
         0.5*(c1-b1*c0)*(c2 - b1*c1 + b1*b1*c0 - b2*c0) +
         (-b1*c2 + b1*b1*c1 - b2*c1 + 2.*b1*b2*c0)/3.)*x*x*x
        )*pow(x,c0);

  x = a_3_i/pi;
  c_i = (1. + (c1-b1*c0)*x +
         0.5*((c1-b1*c0)*(c1-b1*c0) + c2 - b1*c1 + b1*b1*c0 - b2*c0)*x*x +
         ((c1-b1*c0)*(c1-b1*c0)*(c1-b1*c0)/6. +
         0.5*(c1-b1*c0)*(c2 - b1*c1 + b1*b1*c0 - b2*c0) +
         (-b1*c2 + b1*b1*c1 - b2*c1 + 2.*b1*b2*c0)/3.)*x*x*x
        )*pow(x,c0);

  R = c_f/c_i;

/*
  cout<<"L = "<<L<<endl;
  cout<<"a_3_f = "<<a_3_f<<endl;
  cout<<"R = "<<R<<endl;
*/
}

//---------------------------------------------------------------
// input is alpha_3 and masses of fermions at Z-scale in MSbar,
// after threshold corrections
//---------------------------------------------------------------

double MZto2GeV(double M_Z, double a_em, double a_3, VectorXd& mYu,
               VectorXd& mYd, VectorXd& mYe,
               double ex_m_b, double ex_m_c, double ex_m_tau, double ex_m_m,
               double ex_m_e, double& Mb, double& Mc)
{
  std::cout<<"MZto2GeVing!"<<std::endl;
  double P;
  P = 1.;

  int nf;
  const double pi = 3.141592653589793;
  double qed_b, qed_c, qed_sd, qed_u, qed_tau, qed_m, qed_e;
  double a_3_f, Rmb, Rmc, Rm2, tmb;

//---------------------------------------------------------------
// QED corrections:

  qed_tau = 1. + a_em*(4. - 3.*std::log(ex_m_tau*ex_m_tau/(M_Z*M_Z)))/(4.*pi);
  qed_m =   1. + a_em*(4. - 3.*std::log(ex_m_m*ex_m_m/(M_Z*M_Z)))/(4.*pi);
  qed_e =   1. + a_em*(4. - 3.*std::log(ex_m_e*ex_m_e/(M_Z*M_Z)))/(4.*pi);

  qed_b  = 1. + a_em*(1./9.)*(4. - 3.*std::log(ex_m_b*ex_m_b/(M_Z*M_Z)))/(4.*pi);
  qed_c  = 1. + a_em*(4./9.)*(4. - 3.*std::log(ex_m_c*ex_m_c/(M_Z*M_Z)))/(4.*pi);
  qed_sd = 1. + a_em*(1./9.)*(4. - 3.*std::log(2.*2./(M_Z*M_Z)))/(4.*pi);
  qed_u  = 1. + a_em*(4./9.)*(4. - 3.*std::log(2.*2./(M_Z*M_Z)))/(4.*pi);

//---------------------------------------------------------------
// QCD corrections:

//----------------------------
// running from M_Z to m_b

  nf = 5;
      
  run_qcd(nf, a_3, a_3_f, M_Z, ex_m_b, Rmb);

  if (Rmb<0.)
  {
    P = P*10.;
    return P;
  }

// HQEtB factor
// from hep-ph/0004189 eq.(13) neglecting std::logs and quarks lighter than mc

  Mb = 1./(1. - 4.*a_3_f/(3.*pi)
                 - (3019./288. + (pi*pi/6.)*(2.+2.*std::log(2.)/3.-(nf-1.)/3.)
                    -71.*(nf-1)/144. - ZETA_3/6. + (4./3.)*del(ex_m_c/ex_m_b)
                   )*a_3_f*a_3_f/(pi*pi));

// decoupling at flavor thresholds
  tmb = 1. + 89.*a_3_f*a_3_f/(432.*pi*pi);

  a_3 = (1. + 11.*a_3_f*a_3_f/(72.*pi*pi)) * a_3_f;

/*
cout<<"Rmb = "<<Rmb<<endl;
cout<<"tmb = "<<tmb<<endl;
cout<<"a_3 at mb = "<<a_3<<endl;
cout<<"HQEtB= "<<Mb<<endl;
*/

//----------------------------
// running from m_b to m_c

  nf = 4;

  run_qcd(nf, a_3, a_3_f, ex_m_b, ex_m_c, Rmc);

  if (Rmc<0.)
  {
    P = P*10.;
    return P;
  }

// HQEtc factor
  Mc = 1./(1. - 4.*a_3_f/(3.*pi)
                 - (3019./288. + (pi*pi/6.)*(2.+2.*std::log(2.)/3.-(nf-1.)/3.)
                    -71.*(nf-1)/144. - ZETA_3/6.
                   )*a_3_f*a_3_f/(pi*pi));

/*
cout<<"Rmc = "<<Rmc<<endl;
cout<<"a_3 at mc = "<<a_3_f<<endl;
cout<<"HQEtC= "<<Mc<<endl;
*/
           
//----------------------------
// running from m_b to 2 GeV
// for light quarks
 
  nf = 4;

  run_qcd(nf, a_3, a_3_f, ex_m_b, 2., Rm2);

  if (Rm2<0.)
  {
    P = P*10.;
    return P;
  }

/*
cout<<"Rm2 = "<<Rm2<<endl;
cout<<"a_3 at 2 = "<<a_3_f<<endl;
*/

//----------------------------
// corrected fermion masses

  mYd(2) = mYd(2) * Rmb * qed_b;
  mYu(1) = mYu(1) * Rmb * tmb * Rmc * qed_c;

  mYd(1) = mYd(1) * Rmb * tmb * Rm2 * qed_sd;
  mYd(0) = mYd(0) * Rmb * tmb * Rm2 * qed_sd;
  mYu(0) = mYu(0) * Rmb * tmb * Rm2 * qed_u;

  mYe(2) = mYe(2) * qed_tau;
  mYe(1) = mYe(1) * qed_m;
  mYe(0) = mYe(0) * qed_e;

  Mb = mYd(2) * Mb;
  Mc = mYu(1) * Mc;

/*
cout<<"eta_b = "<<Rmb*qed_b<<endl;
cout<<"eta_c = "<<Rmb*tmb * Rmc *qed_c<<endl;
*/

  return P;

}
