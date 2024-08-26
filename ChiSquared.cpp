#include "ChiSquared.h"

// Constructor
ChiSquared::ChiSquared(const std::vector<double>& exp, const std::vector<double>& err)
    : experimental(exp), errors(err) {}

// Operator() to calculate Chi-squared
double ChiSquared::operator()(const std::vector<double>& params) const {
    VectorXcd parametersEW_DR = RunRGEs(params);
    std::vector<double> theoretical = Matching(params, experimental, parametersEW_DR, 246.22);

    double chi2 = 0.0;
    for (size_t i = 0; i < experimental.size(); ++i) {
        double diff = experimental[i] - theoretical[i];
        chi2 += (diff * diff) / (errors[i] * errors[i]);
    }
    return chi2 * RGEsGlobals::penaltyFactor;
}

// Required by Minuit for error definition
double ChiSquared::Up() const {
    return 1.0;
}


