#include "Matching.h"

std::vector<double> Matching(const std::vector<double>& GUTparameters, const std::vector<double>& ExperimentalValues, VectorXcd& parameters, double v, double *chi2)
{
    std::cout<<"Matching!"<<std::endl;
    RGEsGlobals::penaltyFactor = 1;
    double MEW= 91.18;
    std::cout.precision(16);
    MatrixXcd Yu(3,3), Yd(3,3), Ye(3,3),                       // Yukawa matrices
    
      TYu(3,3), TYd(3,3), TYe(3,3),                           // Trilinear soft terms
      MQ2(3,3), ML2(3,3), Mu2(3,3), Md2(3,3), Me2(3,3),        // Soft mass terms (squared)
      mQu(3,3), mQd(3,3),                         // mQ plus D-term contributions (not used)
      Yu_Z(3,3), Yd_Z(3,3), Ye_Z(3,3), Yv_Z(3,3),              // Yukawa matrices @MZ

      TYu_Z(3,3), TYd_Z(3,3), TYe_Z(3,3),                     // Trilinear soft terms @MZ
      MQ2_Z(3,3), ML2_Z(3,3), Mu2_Z(3,3), Md2_Z(3,3), Me2_Z(3,3), // Soft mass terms (squared) @MZ
    
      GuL(6,3), GdL(6,3), GeL(6,3), GvL(3,3),                  // Squark & slepton mixing matrices
      GuR(6,3), GdR(6,3), GeR(6,3), VCKM(3,3), VMNS(3,3),
      U(2,2), V(2,2),                                         // Chargino mixing matrices
      Nt(4,4);                                                // Neutralino mixing matrix

    MatrixXd Nt_mix(4,4);                           // Complex conjugate of Nt
    VectorXcd new_parameters(107);
    VectorXcd M(3);                                      // Gaugino masses

    VectorXd g(3), g_Z(3),                     // Couplings (GUT normalization)
      mYu(3), mYd(3), mYe(3),                 // masses of the fermions
      mMu(6), mMd(6), mMe(6),                 // masses of squarks & sleptons
      mMch(2), mMnt(4),                       // masses of charginos & neutralinos
      mMu_Z(6), mMd_Z(6), mMe_Z(6),
      mMch_Z(2), mMnt_Z(4);

    dcomplex MHu2, MHd2; //Bmu;

    double mHu_new, mHd_new, // calculated from MZ and g
      mstop_1, mstop_2, s2stop, c2stop, sstop, cstop, ssbot, csbot, sstau, cstau,
      msbot_1, msbot_2, s2sbot, c2sbot,
      mstau_1, mstau_2, s2stau, c2stau, Mg, mg_slha; // Gluino mass

    int sign_mu2;
    
    double MW, MZ, sw2, mh;

    double MGUT, tanb, mA, M_Zpole, M_Wpole, mA_pole, chi2_best = 0.0,
      mHp_pole, mH_pole, mh_pole,
      mu2_new, m_A2_new, chi2_EWSB, mA_pole_zal, // mA_pole_zal: used to check convergence of mA_pole, m_A2_new = mA^2
      mHp, mH, c2a, s2a, sw2pole, // mH: CP even heavier Higgs, sw2pole: 1-mW^2/mZ^2
      alpha_em_MZ, vd, vu,
      V_ud, V_us, V_ub, V_cd, V_cs, V_cb, V_td, V_ts, V_tb,
      J, s_del, // J: Jarlskorg, s_del: connected to phase in CKM
      Mb, Mc,   // bottom and charm pole masses
      sw2MSbar; // sw2MSbar to compare with data
    
    double Gmu, alpha_em, alpha3, rho_new;
    double Mt_Pole, Mtau, Ms, Mmu, Mu, Md, Me, Mb_Mc, Md_over_Ms, Q_2, Q;
    
    int k = 0;
    for(int i=0; i<=2; i++)
    {
        for (int j=0; j<=2; j++)
        {
            Yu(i,j) = parameters(k);
            Yd(i,j) = parameters(9+k);
            Ye(i,j) = parameters(18+k);
            TYu(i,j) = parameters(27+k);
            TYd(i,j) = parameters(36+k);
            TYe(i,j) = parameters(45+k);
            MQ2(i,j) = parameters(54+k);
            ML2(i,j) = parameters(63+k);
            Mu2(i,j) = parameters(72+k);
            Md2(i,j) = parameters(81+k);
            Me2(i, j) = parameters(90+k);
            k++;
        }
        M(i) = parameters(99+i);
        g(i) = std::real(parameters(102+i));
    }
    MHu2 = parameters(105);
    MHd2 = parameters(106);
    //Bmu = parameters(107);
    
    Yu = Yu.transpose().eval();     // In model.cpp, Yukawas are defined with doublets on the left,
    Yd = Yd.transpose().eval();     // but RGEs are written with doublets on the right! We Transpose
    Ye = Ye.transpose().eval();     // them here. Furthermore, Yukawas entering RGEs should be complex
                            // conjugate but I think it doesn't make any difference, as far as
    TYu = TYu.transpose().eval();   // we define diagonalization  matrices and VCKM properly, and as far
    TYd = TYd.transpose().eval();   // as soft susy breaking pars are real.
    TYe = TYe.transpose().eval();
    
    tanb = GUTparameters[4];
    
    Yu_Z = Yu;
    Yd_Z = Yd;
    Ye_Z = Ye;
    TYu_Z = TYu;
    TYd_Z = TYd;
    TYe_Z = TYe;
    MQ2_Z = MQ2;
    ML2_Z = ML2;
    Mu2_Z = Mu2;
    Md2_Z = Md2;
    Me2_Z = Me2;
    g_Z = g;
    
    std::vector<double> list = {alpha_em, Gmu, alpha3, MZ, MW, Mu, Md, Mc, Ms, Mt_Pole, Mb, Mtau, Mmu,
        Me, Mb_Mc, Md_over_Ms, Q_2, Q, sw2, V_ud, V_us, V_ub, V_cd, V_cs, V_cb, V_td, V_ts, V_tb, mh};
    //-------------------------------------------------
    // tree level DRbar MAsses, using v = v_best and MA
    //-------------------------------------------------
    double M_Z2 = 0.25 * (g(0) * g(0) * 3./5. + g(1) * g(1)) * v * v;
    MZ  = std::sqrt(M_Z2);
    
    double M_W2 = 0.25 * g(1) * g(1) * v * v;
    MW = std::sqrt(M_W2);
    
    sw2  = g(0) * g(0) * (3./5.) / (g(0) * g(0) * 3./5. + g(1) * g(1));
    
    double cb  = 1. / std::sqrt(1.+tanb*tanb);
    double sb  = tanb * cb;
    
    //Cos(2Beta), Sin(2Beta), Tan(2Beta):
    double c2b = cb * cb - sb * sb;
    double s2b = 2.*sb * cb;
    double t2b = s2b / c2b;

    double mu2 = 0.5 * t2b * (std::real(MHu2) * tanb - std::real(MHd2) / tanb) - M_Z2/2;
    
    if (mu2 < 0.0) {std::cout<<"WARNING! Negative |mu|^2 trial!"<<std::endl; return list;}
    
    double mu = sqrt(fabs(mu2));

    double mA2 = std::real(MHd2) + std::real(MHu2) + 2. * mu2;
    mA2 = std::abs(mA2);
    mA = std::sqrt(mA2);

    double mHp2 = mA2 + M_W2;
    mHp  = sqrt(mHp2);

    double mH2  = 0.5 * (mA2 + M_Z2 + sqrt((mA2+M_Z2) * (mA2+M_Z2) - 4.*mA2 * M_Z2 * c2b * c2b));
    mH   = sqrt(mH2);

    double mh2  = 0.5 * (mA2 + M_Z2 - sqrt((mA2+M_Z2) * (mA2+M_Z2) - 4. * mA2 * M_Z2*c2b * c2b));
    mh = sqrt(mh2);
    
    double t2a   = t2b * (mA2 + M_Z2) / (mA2 - M_Z2);
    c2a   = - c2b * (mA2 - M_Z2) / (mH2 - mh2);
    s2a   = t2a * c2a;
    //To be understood as Tan(2alpha) and so on
    
    mu = GUTparameters[5];
    
    M_Zpole = ExperimentalValues[3];
    M_Wpole = MW;

    mA_pole = mA;
    mA_pole_zal = mA;

    mHp_pole = mHp;
    mH_pole = mH;
    

    //-------------------------------------------------
    // EWSB Loop
    //-------------------------------------------------
    int do_counter = 0;
    do {
      Yu  = Yu_Z;   MQ2 = MQ2_Z;   // Matrices will be altered in
      Yd  = Yd_Z;   Mu2 = Mu2_Z;   // Diogonalize.cpp, so we need
      Ye  = Ye_Z;   Md2 = Md2_Z;   // to restore them before each
      TYu = TYu_Z;  ML2 = ML2_Z;   // run through the loop.
      TYd = TYd_Z;  Me2 = Me2_Z;   // In diag, we also added D terms.
      TYe = TYe_Z;  g   = g_Z;
      double P;
      try{
      P = Diagonalize(v, tanb, mu, g, M, Yu, Yd, Ye, mYu, mYd, mYe, VCKM, TYu, TYd, TYe, MQ2, ML2, Mu2, Md2, Me2, mQu, mQd, mMu, mMd, mMe,
                      GuL, GdL, GeL, GuR, GdR, GeR, GvL, mMch, U, V, mMnt, Nt, Nt_mix, MW, MZ, sw2, mA, mHp, mH, mh, c2a, s2a,
                      mstop_1, mstop_2, s2stop, c2stop, msbot_1, msbot_2, s2sbot, c2sbot, mstau_1, mstau_2, s2stau, c2stau, Mg);

      } catch(...){
        std::cout<<"Too many iterations in Diagonalize: Fixed!"<<std::endl;
        P = 100;
      }
        double epsilon = 1.0e-6;
        
        mMu_Z = mMu;
        mMd_Z = mMd;
        mMe_Z = mMe;
        mMnt_Z = mMnt;
        mMch_Z = mMch;

        
        if (fabs(P-1.)>epsilon) { *chi2 = 1.e13 * P; return list;}

        //-----------------------------------------
        // Integrating out first two family scalars
        //-----------------------------------------
        double invalpha1_un = (4.*M_PI)/(g(0)*g(0));
        double invalpha2_un = (4.*M_PI)/(g(1)*g(1));
        double invalpha3_un = (4.*M_PI)/(g(2)*g(2));

        double invalpha1_c = invalpha1_un - 2./(3.*M_PI)*log(GUTparameters[8]/MEW);
        double invalpha2_c = invalpha2_un - 2./(3.*M_PI)*log(GUTparameters[8]/MEW);
        double invalpha3_c = invalpha3_un - 2./(3.*M_PI)*log(GUTparameters[8]/MEW);

        g(0) = sqrt(4.*M_PI/invalpha1_c);
        g(1) = sqrt(4.*M_PI/invalpha2_c);
        g(2) = sqrt(4.*M_PI/invalpha3_c);

      //Electroweak Symmetry Breaking
      double Pe = EWSB(tanb, mu, mA, std::real(MHu2), std::real(MHd2), g, M, Mg, mYu, mYd, mYe, VCKM, TYu, TYd, TYe, mMu, mMd, mMe, GuL, GdL, GeL, GuR, GdR, GeR, GvL, MQ2, Mu2, Md2, ML2, Me2, mQu, mQd, mstop_1, mstop_2, s2stop, c2stop, msbot_1, msbot_2, s2sbot, c2sbot, mstau_1, mstau_2, s2stau, c2stau, mMch, U, V, mMnt, Nt, mu2_new, m_A2_new, mHu_new, mHd_new, v, MZ, MW, sw2, mHp, mH, mh, c2a, s2a, M_Zpole, M_Wpole, mA_pole, mHp_pole, mH_pole, mh_pole, Gmu, alpha_em, alpha_em_MZ, alpha3, rho_new, sw2pole);

      if (fabs(Pe-1.)>epsilon) { *chi2 = 1.e15 * Pe; return list; }

        sign_mu2 = 1;
        if (mu2_new < 0.) {
          if (m_A2_new < 0.) {
            *chi2 = *chi2 + 1.e3*std::abs(m_A2_new) + 1.e8;
            if (*chi2 < chi2_best) {
              chi2_best = *chi2;
            }
          }
          //seems to work better this way
          sign_mu2 = -1;
          mu2_new = - mu2_new;
        }

        // Check if m_A2_new < 0. mA_pole2 > 0 was tested in ewsb
        if (m_A2_new < 0.) {
          *chi2 = 1.e3*std::abs(m_A2_new) + 1.e8;
          if (*chi2 < chi2_best) {
            chi2_best = *chi2;
          }
          break;
        }

      chi2_EWSB = fabs((mu2_new - mu*mu))/(1000.*1000.) +
        fabs((mA_pole*mA_pole - mA_pole_zal*mA_pole_zal))/
        (1000.*1000.);

      mu = sqrt(mu2_new);
      mA  = sqrt(m_A2_new);
      mA_pole_zal = mA_pole;


    } while (chi2_EWSB > 1.e-8);
    
    //-------------------------------------------------
    // End of EWSB Loop
    //-------------------------------------------------
    
    vu = v / sqrt(2. + 2. /(tanb*tanb));
    vd = v / sqrt(2. + 2.* tanb * tanb);
    
    //Top pole mass:
    Mt_Pole = mYu(2)*( 1. + g(2)*g(2)/(3.*4.*M_PI*M_PI)
      * (6.*std::log(MEW/mYu(2)) + 5.)
      + 10.95 *g(2)*g(2)*g(2)*g(2)/(4.*M_PI*M_PI*4.*M_PI*M_PI)
      + (6.*std::log(MEW/mYu(2)) + 5.) * (4./9.)
      * (g(0)*g(0)*(3./5.)/(g(0)*g(0)*3./5. + g(1)*g(1)))*g(1)*g(1)/(16.*M_PI*M_PI));
    
    V_ud = std::abs(VCKM(0,0));
    V_us = std::abs(VCKM(0,1));
    V_ub = std::abs(VCKM(0,2));
    V_cd = std::abs(VCKM(1,0));
    V_cs = std::abs(VCKM(1,1));
    V_cb = std::abs(VCKM(1,2));
    V_td = std::abs(VCKM(2,0));
    V_ts = std::abs(VCKM(2,1));
    V_tb = std::abs(VCKM(2,2));
    
    double V_ub_over_V_cb = V_ub / V_cb;
    
    /////////////////////////////////////////////////////////////////////////////
    // Running from MZ to 2GeV /////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    // After EWSB, alpha3 == alpha3(MZ)_MSbar while g(3) is not corrected
    // and is still the running DR value at Z scale, which should be used in
    // transforming fermion masses from DR to MS! alpha3 should be compared
    // with exp. value and used for running below the Z scale. alpha_em_MZ is
    // corrected coupling in MS at MZ used for running of fermion masses below MZ
    
    double DRtoMS;
    DRtoMS = 1. + g(2)*g(2)/(3.*4.*M_PI*M_PI);

    mYu(1) = mYu(1) * DRtoMS;
    mYu(0) = mYu(0) * DRtoMS;
    mYd(2) = mYd(2) * DRtoMS;
    mYd(1) = mYd(1) * DRtoMS;
    mYd(0) = mYd(0) * DRtoMS;
    
    // Corrections to lepton masses at weak scale from crossing MZ threshold see
    // Wright eq. (4.4)
    // Should be for all fermions but probably negligible for quarks

    double eta_ewz;

    eta_ewz = 1. + (g(1)*g(1)/(32.*M_PI*M_PI)) *
              ( - ((-3./4.)*(-1.+4.*sw2)*(-1.+4.*sw2) + 1.25)/(4.*(1.-sw2))
                - 0.25 + std::log(MW/MEW)
              );

    mYe(2) = mYe(2) * eta_ewz;
    mYe(1) = mYe(1) * eta_ewz;
    mYe(0) = mYe(0) * eta_ewz;
    
    double P2;
    P2 = MZto2GeV(MEW, alpha_em_MZ, alpha3, mYu, mYd, mYe, ExperimentalValues[10],
                  ExperimentalValues[7], ExperimentalValues[11], ExperimentalValues[12],
                  ExperimentalValues[13], Mb, Mc);

    if (P2!=1.)
    {
      RGEsGlobals::penaltyFactor = 1.e13 * P2;
    }

    Mb = mYd(2);
    Mtau = mYe(2);

    Mc = mYu(1);
    Ms = mYd(1);
    Mmu = mYe(1);

    Mu = mYu(0);
    Md = mYd(0);
    Me = mYe(0);

    Mb_Mc = Mb - Mc;
    Md_over_Ms = Md/Ms;
    Q_2 = (Md-Mu)*(Md+Mu)/( (Ms-0.5*(Md+Mu))*(Ms+0.5*(Md+Mu)) ); // Defined in PDG, better to compare to experiment
    if (Q_2>0.) {
      Q = 1/sqrt(Q_2);
    } else {
      RGEsGlobals::penaltyFactor = 1.e13;
    }
    
    k = 0;
    for(int i=0; i<=2; i++)
    {
        for (int j=0; j<=2; j++)
        {
            parameters(k) = Yu(i,j);
            parameters(9+k) = Yd(i,j);
            parameters(18+k) = Ye(i,j);
            parameters(27+k) = TYu(i,j);
            parameters(36+k) = TYd(i,j);
            parameters(45+k) = TYe(i,j);
            parameters(54+k) = MQ2(i,j);
            parameters(63+k) = ML2(i,j);
            parameters(72+k) = Mu2(i,j);
            parameters(81+k) = Md2(i,j);
            parameters(90+k) = Me2(i, j);
            k++;
        }
        parameters(99+i) = M(i);
        parameters(102+i) = std::real(g(i));
    }
    parameters(105) = MHu2;
    parameters(106) = MHd2;
    //parameters(107) = Bmu;
    
    return list;
}
