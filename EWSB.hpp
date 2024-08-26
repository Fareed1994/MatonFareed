#ifndef EWSB_hpp
#define EWSB_hpp
#include <iostream>
#include <cmath>
#include <complex>
#include "eigen3/Eigen/Dense"

using namespace Eigen;

double EWSB(double tanb, double mu, double mA,
  double mHu, double mHd, VectorXd &g, VectorXcd &M,
  double Mg,
  VectorXd& mYu, VectorXd& mYd, VectorXd& mYe,
  MatrixXcd& VCKM,
  MatrixXcd& TYu, MatrixXcd& TYd,
  MatrixXcd& TYe,
  VectorXd& mMu, VectorXd& mMd, VectorXd& mMe,
  MatrixXcd& GuL, MatrixXcd& GdL,
  MatrixXcd& GeL, MatrixXcd& GuR,
  MatrixXcd& GdR, MatrixXcd& GeR,
  MatrixXcd& GvL,
  MatrixXcd& MQ2, MatrixXcd& Mu2,
  MatrixXcd& Md2, MatrixXcd& ML2,
  MatrixXcd& Me2,
  MatrixXcd& mQu, MatrixXcd& mQd,
  double mstop_1, double mstop_2, double s2stop, double c2stop,
  double msbot_1, double msbot_2, double s2sbot, double c2sbot,
  double mstau_1, double mstau_2, double s2stau, double c2stau,
  VectorXd& mMch,
  MatrixXcd& U, MatrixXcd& V,
  VectorXd& mMnt, MatrixXcd& Nt,
  double& mju2_new, double& m_A2_new, double& mHu_new, double& mHd_new,
  double& v,
  double& M_Z, double& M_W, double& sw2, double& mHp, double& mH,
  double& mh, double& c2a, double& s2a,
  double& M_Zpole, double& M_Wpole, double& mA_pole,
  double& mHp_pole, double& mH_pole, double& mh_pole,
  double& Gmu, double& alpha_em, double& alpha_em_MZ, double& alpha_3,
  double& rho_new, double& sw2pole);

#endif /* EWSB_hpp */
