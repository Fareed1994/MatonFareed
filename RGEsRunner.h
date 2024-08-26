#ifndef RGESRUNNER_H
#define RGESRUNNER_H
#include <iostream>
#include <cmath>
#include "eigen3/Eigen/Dense"

//For linear algebra MatrixXcd, VectorXcd and dcomplex (d stands for double, c for complex):
using namespace Eigen;

namespace RGEsGlobals {
    extern double penaltyFactor;
}

int integrateOdes(VectorXcd &ystart, double from, double to, double eps, double h1, double hmin, VectorXcd (*derivs)(double, const VectorXcd &));
VectorXcd RunRGEs(const std::vector<double>& parameters);

#endif //RGESRUNNER_H
