//--------------------------------------------------------------------------------------------------------------
// Created by Fareed Alasiri on 7/30/24.
//--------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <vector>
#include "eigen3/Eigen/Dense"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnUserParameters.h"
#include "ChiSquared.h"

using namespace ROOT::Minuit2;

int main() {
    std::vector<double> experimental = {137.035999139, 1.1663787e-5, 0.1185, 91.1876, 80.385, (2.3+(0.7-0.5)/2.)*1.e-3, (4.8+(0.5-0.3)/2.)*1.e-3, 1.275, 95e-3, 173.21, 4.18, 1776.86e-3, 105.6583715e-3, 0.510998928e-3, 3.45, 2./(17.+22.), 1./(23.*23.), 23., 0.22333, 0.97425,
        0.2253, 3.85e-3, 0.225, 0.986, 40.8e-3, 0.00840, 0.0400, 1.021, 125.09};
    // Experimental = alpha_em, Gmu, alpha3, MZ, MW, Mu, Md, Mc, Ms, Mt_Pole, Mb, Mtau, Mmu, Me, Mb_Mc, Md_over_Ms, Q_2, Q, sw2, V_ud, V_us, V_ub, V_cd, V_cs, V_cb, V_td, V_ts, V_tb, mh;
    std::vector<double> errors = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    std::vector<double> parameters(12);
    
    parameters[0]  = 25.0;             // 4.*pi / gGUT^2
    parameters[1]  = 2.58374;          // MGUT / 1.e16
    parameters[2]  = 91.0;             // MEW
    parameters[3]  = 0.619038;         // Lambda
    parameters[4]  = 50.0;             // TanBeta
    parameters[5]  = 1055.21;          // mu
    parameters[6]  = 400.0;            // M12
    parameters[7]  = 1.0;              // A
    parameters[8]  = 25000.0;          // M0
    parameters[9]  = 1.93917;          // MHd^2/M0^2
    parameters[10] = 1.60969;          // MHu^2/M0^2
    parameters[11] = -51017.0;         // A0
    
    std::cout.precision(16);
    
    MnUserParameters upar;
    for (size_t i = 0; i < parameters.size(); ++i) {
        upar.Add("Param" + std::to_string(i), parameters[i], 0.1);
    }

    ChiSquared chi2Func(experimental, errors);
    ROOT::Minuit2::MnMigrad migrad(chi2Func, upar);
    FunctionMinimum min = migrad();
    
    if (min.IsValid()) {
        // Access specific properties
        double minValue = min.UserState().Fval();  
        std::cout << "Minimum chi^2: " << minValue << std::endl;
    } else {
        std::cout << "Minimization did not converge or is not valid." << std::endl;
    };
    
    const ROOT::Minuit2::MnUserParameterState& params = min.UserState();
    
    for (unsigned i = 0; i < params.Params().size(); ++i) {
               std::cout << "Parameter " << ": "
                         << params.Value(i) << " +/- " << params.Error(i) << std::endl;
           }
    return 0;
}

