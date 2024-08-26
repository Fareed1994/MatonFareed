#ifndef RGESDEFINER_H
#define RGESDEFINER_H
#include <iostream>
#include <cmath>
#include "eigen3/Eigen/Dense"

//For linear algebra MatrixXcd, VectorXcd and dcomplex (d stands for double, c for complex):
using namespace Eigen;

VectorXcd RGEsystem(double t, const VectorXcd &y);

#endif //RGESDEFINER_H
