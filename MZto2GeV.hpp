#ifndef MZto2GeV_hpp
#define MZto2GeV_hpp
//---------------------------------------------------------
// running of alpha_s and quark masses from MZ to 2 GeV
// using formulas collected in hep-ph/0004189v1
//---------------------------------------------------------

#include <iostream>
#include <cmath>
#include <complex>
#include "eigen3/Eigen/Dense"
#define ZETA_3 1.202057

using namespace Eigen;

double beta0(int nf);
double beta1(int nf);
double beta2(int nf);
double bb1(int nf);
double bb2(int nf);
double gam0(int nf);
double gam1(int nf);
double gam2(int nf);
double del(double x);
void run_qcd(int nf, double a_3_i, double& a_3_f,
             double sc_i, double sc_f, double& R);
double MZto2GeV(double M_Z, double a_em, double a_3, VectorXd& mYu,
               VectorXd& mYd, VectorXd& mYe,
               double ex_m_b, double ex_m_c, double ex_m_tau, double ex_m_m,
               double ex_m_e, double& Mb, double& Mc);

#endif /* MZto2GeV_hpp */
