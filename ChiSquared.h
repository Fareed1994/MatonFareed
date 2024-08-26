#ifndef CHISQUARED_H
#define CHISQUARED_H

#include <vector>
#include "RGEsRunner.h"
#include "Matching.h"
#include "eigen3/Eigen/Dense"
#include "Minuit2/FCNBase.h"

using namespace ROOT::Minuit2;
using Eigen::VectorXd;
using Eigen::VectorXcd;

class ChiSquared : public FCNBase {
public:
    ChiSquared(const std::vector<double>& exp, const std::vector<double>& err);
    virtual ~ChiSquared() = default;

    double operator()(const std::vector<double>& params) const override;
    double Up() const override;

private:
    std::vector<double> experimental;
    std::vector<double> errors;
};

#endif // CHISQUARED_H
