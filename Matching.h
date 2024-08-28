#ifndef MATCHING_h
#include <iostream>
#include <cmath>
#include "ChiSquared.h"
#include "Diagonalization.hpp"
#include "eigen3/Eigen/Dense"
#include "EWSB.hpp"
#include "MZto2GeV.hpp"

using namespace Eigen;

std::vector<double> Matching(const std::vector<double>& GUTparameters, const std::vector<double>& ExperimentalValues, VectorXcd& parameters, double v, double *chi2);

#define MATCHING_h


#endif /* MATCHING_h */
