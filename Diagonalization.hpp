#ifndef Diagonalization_hpp
#define Diagonalization_hpp
#include <iostream>
#include <cmath>
#include <complex>
#include "eigen3/Eigen/Dense"

using namespace Eigen;

double Diagonalize(double v, double tanb, double mu, VectorXd& g, VectorXcd& M,
        MatrixXcd& Yu, MatrixXcd& Yd, MatrixXcd& Ye,
        VectorXd& mYu, VectorXd& mYd, VectorXd& mYe, MatrixXcd& VCKM,
        MatrixXcd& TYu, MatrixXcd& TYd, MatrixXcd& TYe,
        MatrixXcd& MQ2, MatrixXcd& Mu2, MatrixXcd& Md2,
        MatrixXcd& ML2, MatrixXcd& Me2,
        MatrixXcd& mQu, MatrixXcd& mQd,
        VectorXd& mMu, VectorXd& mMd, VectorXd& mMe,
        MatrixXcd& GuL, MatrixXcd& GdL, MatrixXcd& GeL,
        MatrixXcd& GuR, MatrixXcd& GdR, MatrixXcd& GeR,
        MatrixXcd& GvL,
        VectorXd& mMch, MatrixXcd& U, MatrixXcd& V,
        VectorXd& mMnt, MatrixXcd&Nt, MatrixXd& Nt_mix,
        double M_Z, double M_W, double sw2, double mA, double mHp, double mH,
        double mh, double c2a, double s2a,
        double& mstop_1, double& mstop_2, double& s2stop, double& c2stop,
        double& msbot_1, double& msbot_2, double& s2sbot, double& c2sbot,
        double& mstau_1, double& mstau_2, double& s2stau, double& c2stau, double& Mg);


#endif /* Diagonalization_hpp */
