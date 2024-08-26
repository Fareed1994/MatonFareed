//-----------------------------------------------------------------------------
// Calculating threshold corrections to quarks and leptons at Z scale.
//-----------------------------------------------------------------------------
#ifndef ThresholdCorrectionsFermions_hpp
#define ThresholdCorrectionsFermions_hpp
#include <iostream>
#include <cmath>
#include <complex>
#include "eigen3/Eigen/Dense"

using namespace Eigen;

void ThresholdCorrectionsFermions(double v, double tanb, VectorXd &g,
        VectorXcd &M,
        VectorXd& mYu, VectorXd& mYd, VectorXd& mYe,
        MatrixXcd& VCKM,
        MatrixXcd& Yu, MatrixXcd& Yd,
        MatrixXcd& Ye,
        VectorXd& mMu, VectorXd& mMd, VectorXd& mMe,
        MatrixXcd& GuL, MatrixXcd& GdL,
        MatrixXcd& GeL,
        MatrixXcd& GuR, MatrixXcd& GdR,
        MatrixXcd& GeR, MatrixXcd& GvL,
        VectorXd& mMch, MatrixXcd& U,
        MatrixXcd& V,
        VectorXd& mMnt, MatrixXcd& Nt,
        double M_Z, double M_W, double sw2, double mA, double mHp, double mH, double mh, double c2a, double s2a);

#endif /* ThresholdCorrectionsFermions_hpp */
