#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H
#include <iostream>
#include <cmath>
#include "eigen3/Eigen/Dense"


//For linear algebra MatrixXcd, VectorXcd and dcomplex (d stands for double, c for complex):
using namespace Eigen;

int integrateOdes(VectorXcd &ystart, double from, double to, double eps, double h1, double hmin, VectorXcd (*derivs)(double, const VectorXcd &));
int odeStepper(VectorXcd &y, const VectorXcd &dydx, double *x, double htry, double eps, VectorXd &yscal, double *hdid, double *hnext, VectorXd (*derivs)(double, const VectorXcd &));
void rungeKuttaStep(const VectorXcd &y, const VectorXcd &dydx, double x, double h, VectorXcd &yout, VectorXcd &yerr, VectorXcd (*derivs)(double, const VectorXcd &));

#endif //RUNGEKUTTA_H
