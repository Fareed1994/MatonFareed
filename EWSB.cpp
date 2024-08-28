#include "EWSB.hpp"
#include "PassarinoVeltmanFunctions.h"
//------------------------------------------------------------------------------------
// ElectroWeak Symmetry Breaking higgs masses and mu, M_Z and M_W are calculated here
//------------------------------------------------------------------------------------

double rho2(double r) {
  double pi = 3.141592653589793;
  return 19. - 33.*r/2. + 43.*r*r/12. + 7.*r*r*r/120. -
         pi*sqrt(r)* ( 4. - 1.5*r + 3.*r*r/32. + r*r*r/256. ) -
         pi*pi* ( 2. - 2.*r + 0.5*r*r ) - log(r) * (3.*r - 0.5*r*r);
}

double EWSB(double tanb, double mu, double mA,
  double mHu, double mHd, VectorXd &g, VectorXcd &M,
  double Mg,
  VectorXd& mYu, VectorXd& mYd, VectorXd& mYe,
  MatrixXcd& VCKM,
  MatrixXcd& TYu, MatrixXcd& TYd,
  MatrixXcd& TYe,
  VectorXd& mMu, VectorXd& mMd, VectorXd& mMe,
  MatrixXcd& GuL, MatrixXcd& GdL,
  MatrixXcd& GeL, MatrixXcd& GuR,
  MatrixXcd& GdR, MatrixXcd& GeR,
  MatrixXcd& GvL,
  MatrixXcd& MQ2, MatrixXcd& Mu2,
  MatrixXcd& Md2, MatrixXcd& ML2,
  MatrixXcd& Me2,
  MatrixXcd& mQu, MatrixXcd& mQd,
  double mstop_1, double mstop_2, double s2stop, double c2stop,
  double msbot_1, double msbot_2, double s2sbot, double c2sbot,
  double mstau_1, double mstau_2, double s2stau, double c2stau,
  VectorXd& mMch,
  MatrixXcd& U, MatrixXcd& V,
  VectorXd& mMnt, MatrixXcd& Nt,
  double& mju2_new, double& m_A2_new, double& mHu_new, double& mHd_new,
  double& v,
  double& M_Z, double& M_W, double& sw2, double& mHp, double& mH,
  double& mh, double& c2a, double& s2a,
  double& M_Zpole, double& M_Wpole, double& mA_pole,
  double& mHp_pole, double& mH_pole, double& mh_pole,
  double& Gmu, double& alpha_em, double& alpha_em_MZ, double& alpha_3,
  double& rho_new, double& sw2pole)
{
  std::cout<<"EWSBing!"<<std::endl;
  double Pe;
  Pe = 1.;

  const double pi = M_PI;
  const double sixteenpi2 = 16.* pi * pi;
  double MEW= 91.18;

  double M_Z2, M_W2, M_Zpole2, M_Wpole2, v2, vu, vd, cw2, sw, cw;
  double mA2, mHp2, mH2, mh2;
  double p, X;
  double PI_ZZ, PI_ZZ0, PI_ZZ0_new, PI_WW, PI_WW0, PI_WW0_new,
         PI_AA, PI_HH, PI0, PI_junk;
  double PI_11, PI_12, PI_22, PI_11h, PI_12h, PI_22h;
  double mHu2bar, mHd2bar, t1_v1, t2_v2, bA;
  double cb, sb, cb2, sb2, t2b, c2b, s2b;
  double sa2, ca2, t2a, cb_a2, sb_a2, sa, ca, capb, sapb;

  // both sin beta and cos beta are defined positive
  cb  = 1./sqrt(1.+tanb*tanb);
  sb  = tanb*cb;
  cb2 = cb*cb;
  sb2 = 1. - cb2;

  c2b = cb*cb - sb*sb;
  s2b = 2.*sb*cb;
  t2b = s2b/c2b;

  MatrixXcd POM(6,6), POMe(6,3);
  MatrixXcd A(2,2), C(2,2), AB(2,4), CD(2,4);
  double sstop2, cstop2, ssbot2, csbot2, sstau2, cstau2;
  double sstop, cstop, ssbot, csbot, sstau, cstau;
  MatrixXd G_stop(2,2), G_sbot(2,2), G_stau(2,2);

  double sign;
  if (s2stop < 0.) { sign = -1. ;}
   else {sign = 1.;}
  sstop2 = 0.5 - 0.5*c2stop;
  cstop2 = 1. - sstop2;
  sstop = sqrt(sstop2);//*sstop2/fabs(sstop2);
  cstop = sign*sqrt(cstop2);
  G_stop(0,0) = cstop;
  G_stop(0,1) = sstop;
  G_stop(1,0) = -sstop;
  G_stop(1,1) = cstop;

  if (s2sbot < 0.) { sign = -1. ;}
   else {sign = 1.;}
  ssbot2 = 0.5 - 0.5*c2sbot;
  csbot2 = 1. - ssbot2;
  ssbot = sqrt(ssbot2);//*ssbot2/fabs(ssbot2);
  csbot = sign*sqrt(csbot2);
  G_sbot(0,0) = csbot;
  G_sbot(0,1) = ssbot;
  G_sbot(1,0) = -ssbot;
  G_sbot(1,1) = csbot;

  if (s2stau < 0.) { sign = -1. ;}
   else {sign = 1.;}
  sstau2 = 0.5 - 0.5*c2stau;
  cstau2 = 1. - sstau2;
  sstau = sqrt(sstau2);//*sstau2/fabs(sstau2);
  cstau = sign*sqrt(cstau2);
  G_stau(0,0) = cstau;
  G_stau(0,1) = sstau;
  G_stau(1,0) = -sstau;
  G_stau(1,1) = cstau;


  /////////////////////////////////////////////////////////////////////////////
  // tree level DRbar masses, using v = v_best, mu and mA
  /////////////////////////////////////////////////////////////////////////////

  //M_Z2 = 0.25*(g(0)*g(0)*3./5. + g(1)*g(1))*v2;
  //M_Z  = sqrt(M_Z2);
  M_Z2 = M_Z * M_Z;
  //M_W2 = 0.25*g(1)*g(1)*v2;
  //M_W = sqrt(M_W2);
  M_W2 = M_W * M_W;
  sw2 = g(0)*g(0)*(3./5.)/(g(0)*g(0)*3./5. + g(1)*g(1));
  cw2  = 1. - sw2;
  cw  = sqrt(cw2);
  sw  = sqrt(sw2);
  mA2  = mA*mA;
  mHp2 = mHp * mHp;
  mH2  = mH * mH;
  mh2  = mh * mh;

  t2a   = s2a/c2a;
  sa2 = (1. - c2a)/2.;
  ca2 = 1. - sa2;
  cb_a2 = 0.5*(1.+c2a*c2b+s2a*s2b);
  sb_a2 = 1. - cb_a2;
  sa = -sqrt(sa2);
  ca = sqrt(ca2);
  sapb = sa*cb + ca*sb;
  capb = ca*cb - sa*sb;

  /////////////////////////////////////////////////////////////////////////////
  //M_Z self-energy : PIERCE //////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  p = M_Zpole;

  // gauge boson and higgs contribution
  PI_ZZ =
    sb_a2 * ( tB22(p,mA,mH) + tB22(p,M_Z,mh) + M_Z2 * B0(p,M_Z,mh) ) +
    cb_a2 * ( tB22(p,M_Z,mH) + tB22(p,mA,mh) + M_Z2 * B0(p,M_Z,mH) ) -
    2.*cw2*cw2*(2.*p*p + M_W2 - M_Z2*sw2*sw2/cw2 ) * B0(p,M_W,M_W)   +
    (8.*cw2*cw2 + (cw2-sw2)*(cw2-sw2)) * tB22(p,M_W,M_W)             +
    (cw2-sw2)*(cw2-sw2) * tB22(p,mHp,mHp);
  PI_junk = PI_ZZ;

  // quark and lepton contribution
  for (int i=0; i<3; i++) {
    PI_ZZ = PI_ZZ + 3.*(
      ((0.5-2.*sw2/3.)*(0.5-2.*sw2/3.) + 4.*sw2*sw2/9.) * H(p,mYu(i),mYu(i)) -
      4.*(0.5-2.*sw2/3.)*2.*sw2/3.*mYu(i)*mYu(i) * B0(p,mYu(i),mYu(i)) +

      ((-0.5+sw2/3.)*(-0.5+sw2/3.)+sw2*sw2/9.) * H(p,mYd(i),mYd(i)) +
      4.*(-0.5+sw2/3.)*sw2/3.*mYd(i)*mYd(i) * B0(p,mYd(i),mYd(i))) +

      ((-0.5+sw2)*(-0.5+sw2)+sw2*sw2) * H(p,mYe(i),mYe(i)) +
      4.*(-0.5+sw2)*sw2*mYe(i)*mYe(i) * B0(p,mYe(i),mYe(i)) +

      0.5*0.5 * H(p,0.,0.);
  }
  PI_junk = PI_ZZ;

  // 1st 2nd family squark and slepton contribution
  /*
  for (int i=1; i<=2; i++)
    PI_ZZ = PI_ZZ + 4.*3.* (
            (0.5-2.*sw2/3.)*(0.5-2.*sw2/3.)
                           * tB22(p,sqrt(real(mQu(i,i))),sqrt(real(mQu(i,i)))) +
            4.*sw2*sw2/9.  * tB22(p,sqrt(real(Mu2(i,i))),sqrt(real(Mu2(i,i)))) +

            (-0.5+sw2/3.)*(-0.5+sw2/3.)
                           * tB22(p,sqrt(real(mQd(i,i))),sqrt(real(mQd(i,i)))) +
            sw2*sw2/9.     * tB22(p,sqrt(real(Md2(i,i))),sqrt(real(Md2(i,i)))) )

            + 4.* (
            (-0.5+sw2)*(-0.5+sw2)
                           * tB22(p,sqrt(real(ML2(i,i))),sqrt(real(ML2(i,i)))) +
            sw2*sw2        * tB22(p,sqrt(real(Me2(i,i))),sqrt(real(Me2(i,i)))) +

            0.5*0.5        * tB22(p,mMv(i),mMv(i)) );
  */

  // 3rd family squark and slepton contribution
  PI_ZZ = PI_ZZ + 4.*3.* (
    ((0.5-2.*sw2/3.)*cstop2 - 2.*sw2/3.*sstop2)*
    ((0.5-2.*sw2/3.)*cstop2 - 2.*sw2/3.*sstop2)*tB22(p,mstop_1,mstop_1) +
    2.*0.5*0.5*cstop2*sstop2*tB22(p,mstop_1,mstop_2) +
    (2.*sw2/3.*cstop2 - (0.5-2.*sw2/3.)*sstop2)*
    (2.*sw2/3.*cstop2 - (0.5-2.*sw2/3.)*sstop2)*tB22(p,mstop_2,mstop_2) +

    ((-0.5+sw2/3.)*csbot2 + sw2/3.*ssbot2)*
    ((-0.5+sw2/3.)*csbot2 + sw2/3.*ssbot2)*tB22(p,msbot_1,msbot_1) +
    2.*0.5*0.5*csbot2*ssbot2*tB22(p,msbot_1,msbot_2) +
    (-sw2/3.*csbot2 - (-0.5+sw2/3.)*ssbot2)*
    (-sw2/3.*csbot2 - (-0.5+sw2/3.)*ssbot2)*tB22(p,msbot_2,msbot_2) )

    + 4.* (
    ((-0.5+sw2)*cstau2 + sw2*sstau2)*
    ((-0.5+sw2)*cstau2 + sw2*sstau2)*tB22(p,mstau_1,mstau_1) +
    2.*0.5*0.5*cstau2*sstau2*tB22(p,mstau_1,mstau_2) +
    (-sw2*cstau2 - (-0.5+sw2)*sstau2)*
    (-sw2*cstau2 - (-0.5+sw2)*sstau2)*tB22(p,mstau_2,mstau_2) );
  PI_junk = PI_ZZ;

  // neutralino contribution
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_ZZ = PI_ZZ + real( (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))*
                            (Nt(i,2)*conj(Nt(j,2)) - Nt(i,3)*conj(Nt(j,3))) )
                            * H(p,mMnt(i),mMnt(j))/4.
              - 0.5* real( (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))*
                           (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))  )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j)) ;
   PI_junk = PI_ZZ;

  // chargino contribution
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_ZZ = PI_ZZ + (
          real( (cw2*conj(V(i,0))*V(j,0) + 0.5*(cw2-sw2)*conj(V(i,1))*V(j,1))*
                (cw2*V(i,0)*conj(V(j,0)) + 0.5*(cw2-sw2)*V(i,1)*conj(V(j,1))) )
        + real( (cw2*conj(U(i,0))*U(j,0) + 0.5*(cw2-sw2)*conj(U(i,1))*U(j,1))*
                (cw2*U(i,0)*conj(U(j,0)) + 0.5*(cw2-sw2)*U(i,1)*conj(U(j,1))) )
                      ) * H(p,mMch(i),mMch(j)) +4.*
          real( (cw2*conj(U(i,0))*U(j,0) + 0.5*(cw2-sw2)*conj(U(i,1))*U(j,1))*
                (cw2*conj(V(i,0))*V(j,0) + 0.5*(cw2-sw2)*conj(V(i,1))*V(j,1)) )
                * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j)) ;

  // overall factor
  PI_ZZ = PI_ZZ * (g(0)*g(0)*3./5. + g(1)*g(1))/sixteenpi2;


  /////////////////////////////////////////////////////////////////////////////
  //M_Z self-energy p=0 : TOMAS ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // Not used anywhere
  // Should remove first two family squarks

  p = 0.;

  // gauge boson and higgs contribution
  PI_ZZ0 = cb_a2 * M_Z2 * B0(p,mH,M_Z)  +  sb_a2 * M_Z2 * B0(p,mh,M_Z)  +
          2.*sw2*sw2 * M_W2 * B0(p,M_W,M_W)  +
          cb_a2 * tB22(p,mH,M_Z)  +  sb_a2 * tB22(p,mh,M_Z)  +
          sb_a2 * tB22(p,mA,mH)   +  cb_a2 * tB22(p,mA,mh)   +
          (1.-2.*sw2)*(1.-2.*sw2) * ( tB22(p,mHp,mHp) + tB22(p,M_W,M_W) ) +
          (1.- sw2)*(1.- sw2) * ( -4.*A0(M_W) + 8.*B22(p,M_W,M_W)
                                  -(4.*p*p + 2.*M_W2)*B0(p,M_W,M_W) );

  PI_ZZ0_new = cb_a2 * M_Z2 * B0(p,mH,M_Z) +
               cb_a2 * tB22(p,mH,M_Z)  +
               sb_a2 * tB22(p,mA,mH)   +  cb_a2 * tB22(p,mA,mh);

  // quark and lepton contribution
  // up-quarks
  X = 4.*(2./3.)*(2./3.)*sw2*sw2 - 4.*0.5*(2./3.)*sw2 +0.5;
  PI_junk = PI_ZZ0;
  for (int i=0; i<3; i++) {
    PI_ZZ0 = PI_ZZ0 +
             3.*( (2./3.)*X*A0(mYu(i)) + (-p*p/3.+2.*mYu(i)*mYu(i))*X/3.
             + (p*p*X/3. + mYu(i)*mYu(i)*X*2./3. - 0.5*mYu(i)*mYu(i))
             * B0(p,mYu(i),mYu(i)) );
    PI_junk = PI_ZZ0;
  }

  PI_junk = PI_ZZ0;

  // down-quarks
  X = 4.*(1./9.)*sw2*sw2 - 4.*0.5*(1./3.)*sw2 +0.5;
  for (int i=0; i<3; i++) {
    PI_ZZ0 = PI_ZZ0 +
             3.*( (2./3.)*X*A0(mYd(i)) + (-p*p/3.+2.*mYd(i)*mYd(i))*X/3.
             + (p*p*X/3. + mYd(i)*mYd(i)*X*2./3. - 0.5*mYd(i)*mYd(i))
             * B0(p,mYd(i),mYd(i)) );
    PI_junk = PI_ZZ0;
  }
  PI_junk = PI_ZZ0;

  // charged leptons
  X = 4.*sw2*sw2 - 4.*0.5*sw2 +0.5;
  for (int i=0; i<3; i++) {
    PI_ZZ0 = PI_ZZ0 + 1.*( (2./3.)*X*A0(mYe(i)) + (-p*p/3.+2.*mYe(i)*mYe(i))*X/3.
                        + (p*p*X/3. + mYe(i)*mYe(i)*X*2./3. - 0.5*mYe(i)*mYe(i))
                        * B0(p,mYe(i),mYe(i)) );
    PI_junk = PI_ZZ0;
  }
  PI_junk = PI_ZZ0;

  // neutrinos; masses = 0.
  X = 0.5;
  PI_ZZ0 = PI_ZZ0 - p*p*X/3. + p*p*X*B0(p,0.,0.);
  PI_junk = PI_ZZ0;

  // squark and slepton contribution

  PI0 = PI_ZZ0;

  // up-squarks
  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
      POM(a,b) = std::complex<double> (0.,0.);
      for (int i=0; i<3; i++) {
        POM(a,b) = POM(a,b) - GuL(a,i) * conj(GuL(b,i));
      }
    }
    POM(a,a) =  2.*(2./3.)*sw2  + POM(a,a);
  }

  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
      PI_ZZ0 = PI_ZZ0 + 3. * tB22(p, mMu(a), mMu(b)) * real(POM(a,b)*POM(b,a));
    }
  }
  PI_junk = PI_ZZ0;

  // down-squarks
  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
      POM(a,b) = std::complex<double> (0.,0.);
      for (int i=0; i<3; i++) {
        POM(a,b) = POM(a,b) - GdL(a,i) * conj(GdL(b,i));
      }
    }
    POM(a,a) =  2./3.*sw2  + POM(a,a);
  }

  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
      PI_ZZ0 = PI_ZZ0 +  3.* tB22(p, mMd(a), mMd(b)) * std::real(POM(a,b) * POM(b,a));
    }
  }
  PI_junk = PI_ZZ0;

  // charged sleptons
  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
      POM(a,b) = std::complex<double> (0.,0.);
      for (int i=0; i<3; i++) {
        POM(a,b) = POM(a,b) - GeL(a,i) * conj(GeL(b,i));
      }
    }
    POM(a,a) =  2.*sw2  + POM(a,a);
  }

  for (int a=0; a<6; a++) {
    for (int b=0; b<6; b++) {
        PI_ZZ0 = PI_ZZ0 +  tB22(p, mMe(a), mMe(b)) * std::real(POM(a,b) * POM(b,a));
    }
  }
  PI_junk = PI_ZZ0;

  PI_ZZ0_new = PI_ZZ0_new + PI_ZZ0 - PI0;
  PI0 = PI_ZZ0;

  // neutralino contribution
  for (int i=0; i<4; i++) {
    for (int j=0; j<4; j++) {
      PI_ZZ0 = PI_ZZ0 + real( (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))*
                            (Nt(i,2)*conj(Nt(j,2)) - Nt(i,3)*conj(Nt(j,3))) )
                            * H(p,mMnt(i),mMnt(j))/4.
              - 0.5* real( (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))*
                           (conj(Nt(i,2))*Nt(j,2) - conj(Nt(i,3))*Nt(j,3))  )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j)) ;

      PI_ZZ0_new = PI_ZZ0_new + PI_ZZ0 - PI0;
      PI0 = PI_ZZ0;
    }
  }
  PI_junk = PI_ZZ0;

  //Chargino as in Pierce:
  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_ZZ0 = PI_ZZ0 + (
          std::real( (cw2*conj(V(i,0))*V(j,0) + 0.5*(cw2-sw2)*conj(V(i,1))*V(j,1))*
                (cw2*V(i,0)*conj(V(j,0)) + 0.5*(cw2-sw2)*V(i,1)*conj(V(j,1))) )
        + std::real( (cw2*conj(U(i,0))*U(j,0) + 0.5*(cw2-sw2)*conj(U(i,1))*U(j,1))*
                (cw2*U(i,0)*conj(U(j,0)) + 0.5*(cw2-sw2)*U(i,1)*conj(U(j,1))) )
                      ) * H(p,mMch(i),mMch(j)) +4.*
          std::real( (cw2*conj(U(i,0))*U(j,0) + 0.5*(cw2-sw2)*conj(U(i,1))*U(j,1))*
                (cw2*conj(V(i,0))*V(j,0) + 0.5*(cw2-sw2)*conj(V(i,1))*V(j,1)) )
                * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j)) ;

  PI_ZZ0_new = PI_ZZ0_new + PI_ZZ0 - PI0;

  // overall factor
  PI_ZZ0 = PI_ZZ0 * (g(0)*g(0)*3./5. + g(1)*g(1))/sixteenpi2;

  PI_ZZ0_new = PI_ZZ0_new * (g(0)*g(0)*3./5. + g(1)*g(1))/sixteenpi2;

  /////////////////////////////////////////////////////////////////////////////
  //M_W self-energy p = M_Wpole : PIERCE //////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  p = M_Wpole;

  // gauge boson and higgs contribution
  PI_WW = sb_a2 * ( tB22(p,mH,mHp) + tB22(p,mh,M_W) + M_W2 * B0(p,mh,M_W) ) +
          cb_a2 * ( tB22(p,mh,mHp) + tB22(p,mH,M_W) + M_W2 * B0(p,mH,M_W) ) +
          tB22(p,mA,mHp) + (1. + 8.*cw2)* tB22(p,M_Z,M_W) -
          sw2 * ( -8.* tB22(p,M_W,0.) + 4.*p*p* B0(p,M_W,0.) ) -
          ( (4.*p*p + M_Z2 + M_W2)*cw2 - M_Z2*sw2*sw2 ) * B0(p,M_Z,M_W);
  PI_junk = PI_WW;

  // quark and lepton contribution
  for (int i=0; i<3; i++) {
    PI_WW = PI_WW + 0.5*3.*H(p,mYu(i),mYd(i)) + 0.5*H(p,0.,mYe(i));
  }
  PI_junk = PI_WW;

  // 1st 2nd family squark and slepton contribution
  /*
  for (int i=1; i<=2; i++) {
    PI_WW = PI_WW + 2.*3.*tB22(p,sqrt(real(mQu(i,i))),sqrt(real(mQd(i,i))))
                  + 2.*   tB22(p,sqrt(real(ML2(i,i))),mMv(i));
  }
  */

  // 3rd family squark and slepton contribution
  PI_WW = PI_WW + 2.*3.* (
          cstop2*csbot2* tB22(p,mstop_1,msbot_1) +
          cstop2*ssbot2* tB22(p,mstop_1,msbot_2) +
          sstop2*csbot2* tB22(p,mstop_2,msbot_1) +
          sstop2*ssbot2* tB22(p,mstop_2,msbot_2) ) +
          2.* (
          cstau2* tB22(p,mstau_1,0) +
          sstau2* tB22(p,mstau_2,0) );
  PI_junk = PI_WW;

  // chargino - neutralino contribution
  for (int i=0; i<4; i++) {
    for (int j=0; j<2; j++) {
      PI_WW = PI_WW + (
              real( (- conj(Nt(i,1))*V(j,0) + conj(Nt(i,3))*V(j,1)/sqrt(2.))*
                    (- Nt(i,1)*conj(V(j,0)) + Nt(i,3)*conj(V(j,1))/sqrt(2.)) )
            + real( (- conj(Nt(i,1))*U(j,0) - conj(Nt(i,2))*U(j,1)/sqrt(2.))*
                    (- Nt(i,1)*conj(U(j,0)) - Nt(i,2)*conj(U(j,1))/sqrt(2.)) )
                      ) * H(p,mMnt(i),mMch(j)) + 4.*
              real( (- conj(Nt(i,1))*U(j,0) - conj(Nt(i,2))*U(j,1)/sqrt(2.))*
                    (- conj(Nt(i,1))*V(j,0) + conj(Nt(i,3))*V(j,1)/sqrt(2.)) )
              * mMnt(i)*mMch(j) * B0(p,mMnt(i),mMch(j));
    }
  }
  PI_junk = PI_WW;

  // overall factor
  PI_WW = PI_WW * g(1)*g(1)/sixteenpi2;


  /////////////////////////////////////////////////////////////////////////////
  //M_W self-energy p = M_Wpole : PIERCE //////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  p = 0.;

  // gauge boson and higgs contribution
  // later should be p = M_Z - DRbar running ???
  PI_WW0 = sb_a2 * ( tB22(p,mH,mHp) + tB22(p,mh,M_W) + M_W2 * B0(p,mh,M_W) ) +
          cb_a2 * ( tB22(p,mh,mHp) + tB22(p,mH,M_W) + M_W2 * B0(p,mH,M_W) ) +
          tB22(p,mA,mHp) + (1. + 8.*cw2)* tB22(p,M_Z,M_W) -
          sw2 * ( -8.* tB22(p,M_W,0.) + 4.*p*p* B0(p,M_W,0.) ) -
          ( (4.*p*p + M_Z2 + M_W2)*cw2 - M_Z2*sw2*sw2 ) * B0(p,M_Z,M_W);
  PI_junk = PI_WW0;

  // quark and lepton contribution
  for (int i=0; i<3; i++)
    PI_WW0 = PI_WW0 + 0.5*3.*H(p,mYu(i),mYd(i)) + 0.5*H(p,0.,mYe(i));
  PI_junk = PI_WW0;

  // squark and slepton contribution
  PI0 = PI_WW0;

  for (int i=0; i<2; i++)
    PI_WW0 = PI_WW0 + 2.*3.*tB22(p,std::sqrt(std::real(mQu(i,i))),std::sqrt(std::real(mQd(i,i))))
                  + 2.*   tB22(p,std::sqrt(std::real(ML2(i,i))),0);

  PI_WW0 = PI_WW0 + 2.*3.* (
          cstop2*csbot2* tB22(p,mstop_1,msbot_1) +
          cstop2*ssbot2* tB22(p,mstop_1,msbot_2) +
          sstop2*csbot2* tB22(p,mstop_2,msbot_1) +
          sstop2*ssbot2* tB22(p,mstop_2,msbot_2) ) +
          2.* (
          cstau2* tB22(p,mstau_1,0) +
          sstau2* tB22(p,mstau_2,0) );
  PI_junk = PI_WW0;

  PI_WW0_new = PI_WW0_new + PI_WW0 - PI0;
  PI0 = PI_WW0;

//*****************************************************************************8
// chargino - neutralino contribution

  for (int i=0; i<4; i++)
    for (int j=0; j<2; j++)
      PI_WW0 = PI_WW0 + (
              real( (- conj(Nt(i,1))*V(j,0) + conj(Nt(i,3))*V(j,1)/sqrt(2.))*
                    (- Nt(i,1)*conj(V(j,0)) + Nt(i,3)*conj(V(j,1))/sqrt(2.)) )
            + real( (- conj(Nt(i,1))*U(j,0) - conj(Nt(i,2))*U(j,1)/sqrt(2.))*
                    (- Nt(i,1)*conj(U(j,0)) - Nt(i,2)*conj(U(j,1))/sqrt(2.)) )
                      ) * H(p,mMnt(i),mMch(j)) + 4.*
              real( (- conj(Nt(i,1))*U(j,0) - conj(Nt(i,2))*U(j,1)/sqrt(2.))*
                    (- conj(Nt(i,1))*V(j,0) + conj(Nt(i,3))*V(j,1)/sqrt(2.)) )
              * mMnt(i)*mMch(j) * B0(p,mMnt(i),mMch(j));

//  cout<<"PI_WW0 charg-neutr = "<<(PI_WW0-PI_junk) * g(1)*g(1)/sixteenpi2<<"\n";

  PI_WW0_new = PI_WW0_new + PI_WW0 - PI0;

//*****************************************************************************8
// overall factor

  PI_WW0 = PI_WW0 * g(1)*g(1)/sixteenpi2;

//  cout<<"PI_WW0 pierce = "<<PI_WW0<<"\n";

  /////////////////////////////////////////////////////////////////////////////
  // 1 loop MZ and MW pole masses : PIERCE ////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  M_Zpole2 = M_Zpole*M_Zpole;
  M_Z2 = M_Zpole2 + PI_ZZ;

  if (M_Z2 >=0.) M_Z  = sqrt(M_Z2);
  else Pe = Pe*10.;
  if (Pe!=1.) return Pe;

  v2 = (M_Z2) *4./(g(0)*g(0)*3./5. + g(1)*g(1));
  v = sqrt(v2);

  M_W2 = 0.25*g(1)*g(1)*v*v;
  M_W = sqrt(M_W2);
  if (M_W2- PI_WW >=0.) M_Wpole = sqrt(M_W2 - PI_WW);
  else Pe = Pe*10.;
  if (Pe!=1.) return Pe;

  M_Wpole2 = M_Wpole*M_Wpole;

  vu = v / sqrt(2. + 2. /(tanb*tanb));
  vd = v / sqrt(2. + 2.* tanb * tanb);

  /////////////////////////////////////////////////////////////////////////////
  // CP-odd Higgs Self-energy : PIERCE ////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  double PI_AA_pom;

  p = mA_pole;

  PI_AA = 0.;

  for (int i=0; i<3; i++)
    PI_AA = PI_AA +
            3.*cb2* mYu(i)*mYu(i) *
            ( p*p*B0(p,mYu(i),mYu(i)) + 2.*A0(mYu(i)) )/(vu*vu) +
            3.*sb2* mYd(i)*mYd(i) *
            ( p*p*B0(p,mYd(i),mYd(i)) + 2.*A0(mYd(i)) )/(vd*vd) +
               sb2* mYe(i)*mYe(i) *
            ( p*p*B0(p,mYe(i),mYe(i)) + 2.*A0(mYe(i)) )/(vd*vd);
  PI_AA_pom = PI_AA;

  // 3rd family squark and slepton contribution
  PI_AA = PI_AA
              - 3.*( mYu(2)*mYu(2)*cb2/(vu*vu) -
                     g(1)*g(1)*0.5*(0.5 - 2.*sw2/3.)*c2b/cw2 )
              * ( cstop2* A0(mstop_1) + sstop2* A0(mstop_2) )
              - 3.*( mYu(2)*mYu(2)*cb2/(vu*vu) -
                     g(1)*g(1)*0.5*(2.*sw2/3.)*c2b/cw2 )
              * ( sstop2* A0(mstop_1) + cstop2* A0(mstop_2) )

              - 3.*( mYd(2)*mYd(2)*sb2/(vd*vd) -
                     g(1)*g(1)*0.5*(-0.5 + sw2/3.)*c2b/cw2 )
              * ( csbot2* A0(msbot_1) + ssbot2* A0(msbot_2) )
              - 3.*( mYd(2)*mYd(2)*sb2/(vd*vd) +
                     g(1)*g(1)*0.5*(sw2/3.)*c2b/cw2 )
              * ( ssbot2* A0(msbot_1) + csbot2* A0(msbot_2) )

              - ( mYe(2)*mYe(2)*sb2/(vd*vd) -
                     g(1)*g(1)*0.5*(-0.5 + sw2)*c2b/cw2 )
              * ( cstau2* A0(mstau_1) + sstau2* A0(mstau_2) )
              - ( mYe(2)*mYe(2)*sb2/(vd*vd) +
                     g(1)*g(1)*0.5*sw2*c2b/cw2 )
              * ( sstau2* A0(mstau_1) + cstau2* A0(mstau_2) )

              + g(1)*g(1)*0.25*c2b/cw2 * A0(0);

  /*
  // 1st 2nd family squark and slepton contribution
  // Approximation : no mixing in 1st and 2nd family
  for (int i=1; i<=2; i++) {
      PI_AA = PI_AA
              - 3.*( mYu(i)*mYu(i)*cb2/(vu*vu) -
                     g(1)*g(1)*0.5*(0.5 - 2.*sw2/3.)*c2b/cw2 )
              * A0(sqrt(real(mQu(i,i))))
              - 3.*( mYu(i)*mYu(i)*cb2/(vu*vu) -
                     g(1)*g(1)*0.5*(2.*sw2/3.)*c2b/cw2 )
              * A0(sqrt(real(Mu2(i,i))))

              - 3.*( mYd(i)*mYd(i)*sb2/(vd*vd) -
                     g(1)*g(1)*0.5*(-0.5 + sw2/3.)*c2b/cw2 )
              * A0(sqrt(real(mQd(i,i))))
              - 3.*( mYd(i)*mYd(i)*sb2/(vd*vd) +
                     g(1)*g(1)*0.5*(sw2/3.)*c2b/cw2 )
              * A0(sqrt(real(Md2(i,i))))

              - ( mYe(i)*mYe(i)*sb2/(vd*vd) -
                     g(1)*g(1)*0.5*(-0.5 + sw2)*c2b/cw2 )
              * A0(sqrt(real(ML2(i,i))))
              - ( mYe(i)*mYe(i)*sb2/(vd*vd) +
                     g(1)*g(1)*0.5*sw2*c2b/cw2 )
              * A0(sqrt(real(Me2(i,i))))

              + g(1)*g(1)*0.25*c2b/cw2 * A0(mMv(i));
  }
  */

  MatrixXd L_ALR(2,2);

  // Up Squark ////////////////////////////////////////////////////////////////
  L_ALR(0,0) = 0;
  L_ALR(0,1) = ( mYu(2)*mu*sb/vu + std::real(TYu(2,2))*cb )/sqrt(2.);
  L_ALR(1,0) = - L_ALR(0,1);
  L_ALR(1,1) = 0;

  // matrix in u1u2 basis
  // gives the same result as unrotated!!!!
  //L_ALR = G_stop * L_ALR * transpose(G_stop);

  // 3rd family up squarks
  PI_AA = PI_AA + 3.* (
              L_ALR(0,0) * L_ALR(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_ALR(0,1) * L_ALR(0,1) * B0(p,mstop_1,mstop_2) +
              L_ALR(1,1) * L_ALR(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // 1st 2nd family up squarks
  for (int i=1; i<=2; i++) {
    L_ALR(1,1) = 0;
    L_ALR(1,2) = ( mYu(i)*mu*sb/vu + real(TYu(i,i))*cb )/sqrt(2.);
    L_ALR(2,1) = - L_ALR(1,2);
    L_ALR(2,2) = 0;

    PI_AA = PI_AA + 3.* 2.*  L_ALR(1,2) * L_ALR(1,2)
                    * B0(p,sqrt(real(mQu(i,i))),sqrt(real(Mu2(i,i))));
  }
  */

  // Down Squark //////////////////////////////////////////////////////////////
  L_ALR(0,0) = 0;
  L_ALR(0,1) = ( mYd(2)*mu*cb/vd + std::real(TYd(2,2))*sb )/sqrt(2.);
  L_ALR(1,0) = - L_ALR(0,1);
  L_ALR(1,1) = 0;

  // matrix in u1u2 basis
  //L_ALR = G_sbot * L_ALR * transpose(G_sbot);

  // 3rd family down squarks
  PI_AA = PI_AA + 3.* (
              L_ALR(0,0) * L_ALR(0,0) * B0(p,msbot_1,msbot_1) +
          2.* L_ALR(0,1) * L_ALR(0,1) * B0(p,msbot_1,msbot_2) +
              L_ALR(1,1) * L_ALR(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // 1st 2nd family down squarks
  for (int i=1; i<=2; i++) {
    L_ALR(1,1) = 0;
    L_ALR(1,2) = ( mYd(i)*mu*cb/vd + real(TYd(i,i))*sb )/sqrt(2.);
    L_ALR(2,1) = - L_ALR(1,2);
    L_ALR(2,2) = 0;

    PI_AA = PI_AA + 3.* 2.*  L_ALR(1,2) * L_ALR(1,2)
                    * B0(p,sqrt(real(mQd(i,i))),sqrt(real(Md2(i,i))));
  }
  */

  // Charged Sleptons /////////////////////////////////////////////////////////
  L_ALR(0,0) = 0;
  L_ALR(0,1) = ( mYe(2)*mu*cb/vd + real(TYe(2,2))*sb )/sqrt(2.);
  L_ALR(1,0) = - L_ALR(0,1);
  L_ALR(1,1) = 0;

  // matrix in u1u2 basis
  //L_ALR = G_stau * L_ALR * transpose(G_stau);

    PI_AA = PI_AA +
                L_ALR(0,0) * L_ALR(0,0) * B0(p,mstau_1,mstau_1) +
            2.* L_ALR(0,1) * L_ALR(0,1) * B0(p,mstau_1,mstau_2) +
                L_ALR(1,1) * L_ALR(1,1) * B0(p,mstau_2,mstau_2) ;

  /*
  for (int i=1; i<=2; i++) {
    L_ALR(1,1) = 0;
    L_ALR(1,2) = ( mYe(i)*mu*cb/vd + real(TYe(i,i))*sb )/sqrt(2.);
    L_ALR(2,1) = - L_ALR(1,2);
    L_ALR(2,2) = 0;

    //PI_AA = PI_AA + 2.*  L_ALR(1,2) * L_ALR(1,2)
    //                * B0(p,sqrt(real(ML2(i,i))),sqrt(real(Me2(i,i))));
  }
  */

//cout<<"PI_AA sc = "<<(PI_AA - PI_AA_pom)/(16.*pi*pi)<<endl;
PI_AA_pom = PI_AA;

// other Higgses

  PI_AA = PI_AA + 0.25 * g(1)*g(1) * (
          2.* F(p,mHp,M_W) + sb_a2* F(p,mH,M_Z)/cw2
          + cb_a2* F(p,mh,M_Z)/cw2 );


  PI_AA = PI_AA + 0.25 * g(1)*g(1)/cw2 * (
          c2b*capb*c2b*capb * M_Z2 * B0(p,mA,mH)  +
          c2b*sapb*c2b*sapb * M_Z2 * B0(p,mA,mh)  +
          s2b*capb*s2b*capb * M_Z2 * B0(p,M_Z,mH) +
          s2b*sapb*s2b*sapb * M_Z2 * B0(p,M_Z,mh)
          - 0.5 * 3.*c2b*c2b * A0(mA) - 0.5* (3.*s2b*s2b - 1.)
          * A0(M_Z)
          + 0.5 *c2b*c2a * A0(mH) - 0.5 *c2b*c2a * A0(mh) );

  PI_AA = PI_AA + 0.5* g(1)*g(1) * M_W2 * B0(p,M_W,mHp)
          - 0.25 * g(1)*g(1)/cw2 * ( c2b*c2b* A0(mHp) +
          (cw2*(1. + s2b*s2b) - sw2*c2b*c2b) * A0(M_W) )
          - 2.*g(1)*g(1)* A0(M_W) - g(1)*g(1)* A0(M_Z)/cw2;



// neutralino contribution

  MatrixXcd AnnA(4,4), BnnA(4,4);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      AnnA(i,j) = std::complex<double> (0.,0.);

  AnnA(0,2) = -0.5* g(0)*sqrt(3./5.) * (-sb) ;
  AnnA(2,0) = AnnA(0,2);
  AnnA(0,3) = -0.5* g(0)*sqrt(3./5.) * cb;
  AnnA(3,0) = AnnA(0,3);

  AnnA(1,2) = 0.5* g(1) * (-sb);
  AnnA(2,1) = AnnA(1,2);
  AnnA(1,3) = 0.5* g(1) * cb;
  AnnA(3,1) = AnnA(1,3);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      BnnA(i,j) = - AnnA(i,j);

// to mass eigenstates basis

  AnnA = (Nt.conjugate() * AnnA * Nt.adjoint()).eval();

  BnnA = (Nt * BnnA * Nt.transpose()).eval();

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_AA = PI_AA + 0.5*
              real( AnnA(i,j)*conj(AnnA(i,j)) + BnnA(i,j)*conj(BnnA(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              2.* real( AnnA(i,j)*conj(BnnA(i,j)) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

// chargino contribution

  MatrixXcd AccA(2,2), BccA(2,2);

  AccA(0,0) = std::complex<double> (0.,0.);
  AccA(1,1) = std::complex<double> (0.,0.);

  AccA(0,1) = g(1) * (-sb) / sqrt(2.) ;

  AccA(1,0) = - g(1) * cb / sqrt(2.) ;

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      BccA(i,j) = - AccA(j,i);

// to mass eigenstates basis

  AccA = (V.conjugate() * AccA * U.adjoint()).eval();

  BccA = (U * BccA * V.transpose()).eval();

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_AA = PI_AA +
              real( AccA(i,j)*conj(AccA(i,j)) + BccA(i,j)*conj(BccA(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              4.* real( AccA(i,j)*conj(BccA(i,j)) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//*****************************************************************************8
// overall factor

   PI_AA = PI_AA/(16.*pi*pi);

//cout<<"PI_AA = "<<PI_AA<<endl;

  /////////////////////////////////////////////////////////////////////////////
  // Charged Higgs self-energy // Added by ZJ /////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  PI_HH = 0.;
  /*
  p = mHp_pole;

  PI_HH = 0.;

  // quark contribution
  for (int i=1; i<=3; i++) {
    PI_HH = PI_HH
      + 3 * ((mYu(i)*mYu(i)/(vu*vu)*cb2 + mYd(i)*mYd(i)/(vd*vd)*sb2)
      * G(p, mYu(i), mYd(i)) - 2*mYu(i)/vu*mYd(i)/vd*mYu(i)*mYd(i)*s2b
      * B0(p, mYu(i), mYd(i)));
  }
  // lepton contribution
  for (int i=1; i<=3; i++) {
    PI_HH = PI_HH +
      (0 + mYe(i)*mYe(i)/(vd*vd)*sb2) * G(p, 0, mYe(i)) - 0;
  }

  // 1st & 2nd generation squark contribution
  // Do not rotate to 12 basis: Mixing is small for 1st & 2nd generation
  MatrixXd L_HLRq(2,2);
  for (int i=1; i<=2; i++) {
    L_HLRq(1,1) = g(1)*M_W*s2b/sqrt(2) - mYu(i)/vu*mYu(i)*cb -
      mYd(i)/vd*mYd(i)*sb;
    L_HLRq(1,2) = mYd(i)/vd*mu*cb - real(TYd(i,i))*sb;
    L_HLRq(2,1) = mYu(i)/vu*mu*sb - real(TYu(i,i))*cb;
    L_HLRq(2,2) = -mYu(i)/vu*mYd(i)*cb - mYd(i)/vd*mYu(i)*sb;
    PI_HH = PI_HH
      + 3 * L_HLRq(1,1)*L_HLRq(1,1) * B0(p,sqrt(real(mQu(i,i))),
        sqrt(real(mQd(i,i))))
      + 3 * L_HLRq(1,2)*L_HLRq(1,2) * B0(p,sqrt(real(mQu(i,i))),
        sqrt(real(Md2(i,i))))
      + 3 * L_HLRq(2,1)*L_HLRq(2,1) * B0(p,sqrt(real(Mu2(i,i))),
        sqrt(real(mQd(i,i))))
      + 3 * L_HLRq(2,2)*L_HLRq(2,2) * B0(p,sqrt(real(Mu2(i,i))),
        sqrt(real(Md2(i,i))));
  }
  // 1st & 2nd generation slepton contribution
  MatrixXd L_HLRl(2,2);
  for (int i=1; i<=2; i++) {
    L_HLRl(1,1) = g(1)*M_W*s2b/sqrt(2) - 0 - mYe(i)/vd*mYe(i)*sb;
    L_HLRl(1,2) = mYe(i)/vd*mu*cb - real(TYe(i,i))*sb;
    L_HLRl(2,1) = 0;
    L_HLRl(2,2) = 0;
    PI_HH = PI_HH
      + L_HLRl(1,1)*L_HLRl(1,1) * B0(p,0,sqrt(real(ML2(i,i))))
      + L_HLRl(1,2)*L_HLRl(1,2) * B0(p,0,sqrt(real(Me2(i,i))))
      + L_HLRl(2,1)*L_HLRl(2,1) * B0(p,0,sqrt(real(ML2(i,i))))
      + L_HLRl(2,2)*L_HLRl(2,2) * B0(p,0,sqrt(real(Me2(i,i))));
  }
  // 3rd generation squark contribution
  MatrixXd L_HLRq3(2,2);
  L_HLRq3(1,1) = g(1)*M_W*s2b/sqrt(2) - mYu(2)/vu*mYu(2)*cb -
    mYd(2)/vd*mYd(2)*sb;
  L_HLRq3(1,2) = mYd(2)/vd*mu*cb - real(TYd(3,3))*sb;
  L_HLRq3(2,1) = mYu(2)/vu*mu*sb - real(TYu(3,3))*cb;
  L_HLRq3(2,2) = -mYu(2)/vu*mYd(2)*cb - mYd(2)/vd*mYu(2)*sb;
  // Rotate to 12 basis
  MatrixXd L_H12q3(2,2);
  L_H12q3 = G_stop * L_HLRq3 * transpose(G_sbot);
  PI_HH = PI_HH
    + 3 * L_H12q3(1,1)*L_H12q3(1,1) * B0(p,mstop_1,msbot_1)
    + 3 * L_H12q3(1,2)*L_H12q3(1,2) * B0(p,mstop_1,msbot_2)
    + 3 * L_H12q3(2,1)*L_H12q3(2,1) * B0(p,mstop_2,msbot_1)
    + 3 * L_H12q3(2,2)*L_H12q3(2,2) * B0(p,mstop_2,msbot_2);
  // 3rd generation slepton contribution
  MatrixXd L_HLRl3(2,2);
  L_HLRl3(1,1) = g(1)*M_W*s2b/sqrt(2) - 0 - mYe(2)/vd*mYe(2)*sb;
  L_HLRl3(1,2) = mYe(2)/vd*mu*cb - real(TYe(3,3))*sb;
  L_HLRl3(2,1) = 0;
  L_HLRl3(2,2) = 0;
  // Rotate to 12 basis
  MatrixXd L_H12l3(2,2);
  MatrixXd NuLRto12(2,2);
  NuLRto12 = 0;
  NuLRto12(1,1) = 1;
  L_H12l3 = NuLRto12 * L_HLRq3 * transpose(G_stau);
  PI_HH = PI_HH
    + 3 * L_H12l3(1,1)*L_H12l3(1,1) * B0(p,0,mstau_1)
    + 3 * L_H12l3(1,2)*L_H12l3(1,2) * B0(p,0,mstau_2)
    + 3 * L_H12l3(2,1)*L_H12l3(2,1) * B0(p,0,mstau_1)
    + 3 * L_H12l3(2,2)*L_H12l3(2,2) * B0(p,0,mstau_2);

  for (int i=1; i<=3; i++) {

  }
  */

//*****************************************************************************8
// CP even Higgs self-energies

// for light Higgs:

  p = mh_pole;

  PI_11 = 0.;

  // 1st line in Pierce ///////////////////////////////////////////////////////
  // Quarks and leptons contribution
  for (int i=0; i<3; i++)
    PI_11 = PI_11 +
            3.*mYd(i)*mYd(i)*( ( p*p-4.*mYd(i)*mYd(i) )* B0(p,mYd(i),mYd(i))
            + 2.* A0(mYd(i)) )/(vd*vd) +

            mYe(i)*mYe(i)*( ( p*p-4.*mYe(i)*mYe(i) )* B0(p,mYe(i),mYe(i))
            + 2.* A0(mYe(i)) )/(vd*vd);

  // 3rd family squarks and sleptons contribution
  PI_11 = PI_11 - 3.*mYd(2)*mYd(2)* ( A0(msbot_1) + A0(msbot_2) )/(vd*vd)
                -    mYe(2)*mYe(2)* ( A0(mstau_1) + A0(mstau_2) )/(vd*vd);

  // 2nd line in Pierce ///////////////////////////////////////////////////////
  // 3rd family squarks and sleptons contribution
  PI_11 = PI_11 - 3.* g(1)*g(1)/(2.*cw2) * (
          (0.5 - 2.*sw2/3.)* ( cstop2* A0(mstop_1) + sstop2* A0(mstop_2) ) +
          2.*sw2/3.        * ( sstop2* A0(mstop_1) + cstop2* A0(mstop_2) ) +

          (-0.5 + sw2/3.)  * ( csbot2* A0(msbot_1) + ssbot2* A0(msbot_2) ) -
          sw2/3.           * ( ssbot2* A0(msbot_1) + csbot2* A0(msbot_2) ) )

                -     g(1)*g(1)/(2.*cw2) * (
          (-0.5 + sw2)     * ( cstau2* A0(mstau_1) + sstau2* A0(mstau_2) ) -
          sw2              * ( sstau2* A0(mstau_1) + cstau2* A0(mstau_2) ) +

          0.5 * A0(0) );

  // 1st 2nd family squarks and sleptons contribution
  /*
  for (int i=1; i<=2; i++)
    PI_11 = PI_11 - 3.* g(1)*g(1)/(2.*cw2) * (
            (0.5 - 2.*sw2/3.)*  A0(sqrt(real(mQu(i,i)))) +
            2.*sw2/3.        *  A0(sqrt(real(Mu2(i,i)))) +

            (-0.5 + sw2/3.)  *  A0(sqrt(real(mQd(i,i)))) -
            sw2/3.           *  A0(sqrt(real(Md2(i,i)))) )

                  -     g(1)*g(1)/(2.*cw2) * (
            (-0.5 + sw2)     *  A0(sqrt(real(ML2(i,i)))) -
            sw2              *  A0(sqrt(real(Me2(i,i)))) +

            0.5 * A0(mMv(i)) );
  */

  // 3rd line in Pierce ///////////////////////////////////////////////////////
  MatrixXd L_1LRu1(2,2), L_1LRu2(2,2), L_1LRu3(2,2),
                 L_1LRd1(2,2), L_1LRd2(2,2), L_1LRd3(2,2),
                 L_1LRe1(2,2), L_1LRe2(2,2), L_1LRe3(2,2);

  /*
  // Up squark : 1st generation
  L_1LRu1(1,1) = g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw;
  L_1LRu1(1,2) = -mYu(0)*mu/(vu*sqrt(2.));
  L_1LRu1(2,1) = L_1LRu1(1,2);
  L_1LRu1(2,2) = g(1)*M_Z*cb*2.*sw2/(3.*cw);

  PI_11 = PI_11 + 3.*(
     L_1LRu1(1,1)*L_1LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
     L_1LRu1(2,2)*L_1LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))));

  // Up squark : 2nd generation
  L_1LRu2(1,1) = g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw;
  L_1LRu2(1,2) = -mYu(1)*mu/(vu*sqrt(2.));
  L_1LRu2(2,1) = L_1LRu2(1,2);
  L_1LRu2(2,2) = g(1)*M_Z*cb*2.*sw2/(3.*cw);

  PI_11 = PI_11 + 3.*(
     L_1LRu2(1,1)*L_1LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
     L_1LRu2(2,2)*L_1LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))));
  */

  // Up squark : 3rd generation
  L_1LRu3(0,0) = g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw;
  L_1LRu3(0,1) = -mYu(2)*mu/(vu*sqrt(2.));
  L_1LRu3(1,0) = L_1LRu3(0,1);
  L_1LRu3(1,1) = g(1)*M_Z*cb*2.*sw2/(3.*cw);

  // matrix in u1u2 basis
  L_1LRu3 = G_stop * L_1LRu3 * G_stop.transpose().eval();

  PI_11 = PI_11 + 3.* (
              L_1LRu3(0,0) * L_1LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_1LRu3(0,1) * L_1LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_1LRu3(1,1) * L_1LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st generation
  L_1LRd1(1,1) = g(1)*M_Z*cb*(-0.5+sw2/3.)/cw + sqrt(2.)*mYd(0)*mYd(0)/vd;
  L_1LRd1(1,2) = real(TYd(1,1))/sqrt(2.);
  L_1LRd1(2,1) = L_1LRd1(1,2);
  L_1LRd1(2,2) = g(1)*M_Z*cb*(-sw2/3.)/cw + sqrt(2.)*mYd(0)*mYd(0)/vd;

  PI_11 = PI_11 + 3.* (
     L_1LRd1(1,1)*L_1LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
     L_1LRd1(2,2)*L_1LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))) );

  // Down squark: 2nd generation
  L_1LRd2(1,1) = g(1)*M_Z*cb*(-0.5+sw2/3.)/cw + sqrt(2.)*mYd(1)*mYd(1)/vd;
  L_1LRd2(1,2) = real(TYd(2,2))/sqrt(2.);
  L_1LRd2(2,1) = L_1LRd2(1,2);
  L_1LRd2(2,2) = g(1)*M_Z*cb*(-sw2/3.)/cw + sqrt(2.)*mYd(1)*mYd(1)/vd;

  PI_11 = PI_11 + 3.* (
          L_1LRd2(1,1)*L_1LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
          L_1LRd2(2,2)*L_1LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))) );
  */

  // Down squark: 3rd generation
  L_1LRd3(0,0) = g(1)*M_Z*cb*(-0.5+sw2/3.)/cw + sqrt(2.)*mYd(2)*mYd(2)/vd;
  L_1LRd3(0,1) = std::real(TYd(2,2))/sqrt(2.);
  L_1LRd3(1,0) = L_1LRd3(0,1);
  L_1LRd3(1,1) = g(1)*M_Z*cb*(-sw2/3.)/cw + sqrt(2.)*mYd(2)*mYd(2)/vd;

  // matrix in u1u2 basis
  L_1LRd3 = G_sbot * L_1LRd3 * G_sbot.transpose().eval();

  PI_11 = PI_11 + 3.* (
              L_1LRd3(0,0) * L_1LRd3(0,0) * B0(p,msbot_1,msbot_1) +
          2.* L_1LRd3(0,1) * L_1LRd3(0,1) * B0(p,msbot_1,msbot_2) +
              L_1LRd3(1,1) * L_1LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Sleptons: 1st generation
  L_1LRe1(1,1) = g(1)*M_Z*cb*(-0.5+sw2)/cw + sqrt(2.)*mYe(0)*mYe(0)/vd;
  L_1LRe1(1,2) = real(TYe(1,1))/sqrt(2.);
  L_1LRe1(2,1) = L_1LRe1(1,2);
  L_1LRe1(2,2) = g(1)*M_Z*cb*(-sw2)/cw + sqrt(2.)*mYe(0)*mYe(0)/vd;

  PI_11 = PI_11 + (
     L_1LRe1(1,1)*L_1LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
     L_1LRe1(2,2)*L_1LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))) );

  // Sleptons: 2nd generation
  L_1LRe2(1,1) = g(1)*M_Z*cb*(-0.5+sw2)/cw + sqrt(2.)*mYe(1)*mYe(1)/vd;
  L_1LRe2(1,2) = real(TYe(2,2))/sqrt(2.);
  L_1LRe2(2,1) = L_1LRe2(1,2);
  L_1LRe2(2,2) = g(1)*M_Z*cb*(-sw2)/cw + sqrt(2.)*mYe(1)*mYe(1)/vd;

  PI_11 = PI_11 + (
     L_1LRe2(1,1)*L_1LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
     L_1LRe2(2,2)*L_1LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))) );
  */

  // Sleptons: 3rd generation
  L_1LRe3(0,0) = g(1)*M_Z*cb*(-0.5+sw2)/cw + sqrt(2.)*mYe(2)*mYe(2)/vd;
  L_1LRe3(0,1) = std::real(TYe(2,2))/sqrt(2.);
  L_1LRe3(1,0) = L_1LRe3(0,1);
  L_1LRe3(1,1) = g(1)*M_Z*cb*(-sw2)/cw + sqrt(2.)*mYe(2)*mYe(2)/vd;

  // matrix in u1u2 basis
  L_1LRe3 = G_stau * L_1LRe3 * G_stau.transpose().eval();

  PI_11 = PI_11 + (
              L_1LRe3(0,0) * L_1LRe3(0,0) * B0(p,mstau_1,mstau_1) +
          2.* L_1LRe3(0,1) * L_1LRe3(0,1) * B0(p,mstau_1,mstau_2) +
              L_1LRe3(1,1) * L_1LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino contribution
  for (int i=0; i<3; i++)
    PI_11 = PI_11 + (g(1)*M_Z*cb*(0.5)/cw)*(g(1)*M_Z*cb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_11 = PI_11 + 0.25 * g(1)*g(1) * (
          sb2* ( 2.* F(p,mHp,M_W) + F(p,mA,M_Z)/cw2 )+
          cb2* ( 2.* F(p,M_W,M_W) + F(p,M_Z,M_Z)/cw2 ) );

  PI_11 = PI_11 + (7./4.)*g(1)*g(1)*cb2* (
          2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 );

  PI_11 = PI_11 - g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2);

//*****************************************************************************8

  PI_11 = PI_11 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (cb*(3.*ca2-sa2) - sb*s2a)*(cb*(3.*ca2-sa2) - sb*s2a)*B0(p,mH,mH) +
          (cb*(3.*sa2-ca2) + sb*s2a)*(cb*(3.*sa2-ca2) + sb*s2a)*B0(p,mh,mh) +
          2.*(-2.*cb*s2a - sb*c2a)*(-2.*cb*s2a - sb*c2a)*B0(p,mH,mh)        +
          c2b*cb*c2b*cb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )                   +
          2.*s2b*cb*s2b*cb * B0(p,M_Z,mA) );

  PI_11 = PI_11 - 0.125*g(1)*g(1)/cw2 * (
          (3.*ca2-sa2) * A0(mH) +
          (3.*sa2-ca2) * A0(mh) +
          c2b * A0(M_Z) - c2b * A0(mA) );


//*****************************************************************************8

  PI_11 = PI_11 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          c2b*cb*c2b*cb * B0(p,M_W,M_W) +
          (-c2b*cb+2.*cw2*cb)*(-c2b*cb+2.*cw2*cb) * B0(p,mHp,mHp) +
          2.*(-s2b*cb+cw2*sb)*(-s2b*cb+cw2*sb) * B0(p,M_W,mHp) );

  PI_11 = PI_11 - 0.25*g(1)*g(1)/cw2 * (
          (cw2 + sw2*c2b) * A0(M_W) +
          (cw2 - sw2*c2b) * A0(mHp) );

//cout<<"PI_11 h= "<<PI_11<<endl;
//*****************************************************************************8
// neutralino contribution

  MatrixXcd Ann1(4,4), Bnn1(4,4);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      Ann1(i,j) = std::complex<double> (0.,0.);
      Bnn1(i,j) = std::complex<double> (0.,0.);
    }

  Ann1(0,2) = -0.5* g(0)*sqrt(3./5.);
  Ann1(2,0) = Ann1(0,2);
  Ann1(1,2) = 0.5* g(1);
  Ann1(2,1) = Ann1(1,2);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      Bnn1(i,j) = Ann1(i,j);

// to mass eigenstates basis

  Ann1 = (Nt.conjugate() * Ann1 * Nt.adjoint()).eval();

  Bnn1 = (Nt * Bnn1 * Nt.transpose()).eval();

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_11 = PI_11 + 0.5*
              real( Ann1(i,j)*conj(Ann1(i,j)) + Bnn1(i,j)*conj(Bnn1(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann1(i,j)*conj(Bnn1(i,j)) +  conj(Ann1(i,j))*Bnn1(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//cout<<"PI_11 n= "<<PI_11<<endl;

//*****************************************************************************8
// chargino contribution

  MatrixXcd Acc1(2,2), Bcc1(2,2);

  Acc1(0,0) = std::complex<double> (0.,0.);
  Acc1(1,1) = std::complex<double> (0.,0.);
  Acc1(0,1) = g(1) / sqrt(2.) ;
  Acc1(1,0) = std::complex<double> (0.,0.);;

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      Bcc1(i,j) =  Acc1(j,i);

// to mass eigenstates basis

  Acc1 = (V.conjugate() * Acc1 * U.adjoint()).eval();

  Bcc1 = (U * Bcc1 * V.transpose()).eval();

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_11 = PI_11 +
              std::real( Acc1(i,j)*conj(Acc1(i,j)) + Bcc1(i,j)*conj(Bcc1(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* std::real( Acc1(i,j)*conj(Bcc1(i,j)) + conj(Acc1(i,j))*Bcc1(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//cout<<"PI_11 c= "<<PI_11<<endl;
//*****************************************************************************8
// overall factor

   PI_11 = PI_11/(16.*pi*pi);

//cout<<"PI_11 = "<<PI_11<<endl;


//*****************************************************************************8
// CP even Higgs self-energies 22

  //p = mh_pole;

  PI_22 = 0.;

  // 1st line in Pierce ///////////////////////////////////////////////////////
  // Quarks and leptons contribution
  for (int i=0; i<3; i++)
    PI_22 = PI_22 +
            3.*mYu(i)*mYu(i)*( ( p*p-4.*mYu(i)*mYu(i) )* B0(p,mYu(i),mYu(i))
            + 2.* A0(mYu(i)) )/(vu*vu) ;

  // 3rd family squarks and sleptons contribution
  PI_22 = PI_22 - 3.*mYu(2)*mYu(2)* ( A0(mstop_1) + A0(mstop_2) )/(vu*vu);

  // 2nd line in Pierce ///////////////////////////////////////////////////////
  // 3rd family squarks and sleptons contribution0
  PI_22 = PI_22 + 3.* g(1)*g(1)/(2.*cw2) * (
          (0.5 - 2.*sw2/3.)* ( cstop2* A0(mstop_1) + sstop2* A0(mstop_2) ) +
          2.*sw2/3.        * ( sstop2* A0(mstop_1) + cstop2* A0(mstop_2) ) +

          (-0.5 + sw2/3.)  * ( csbot2* A0(msbot_1) + ssbot2* A0(msbot_2) ) -
          sw2/3.           * ( ssbot2* A0(msbot_1) + csbot2* A0(msbot_2) ) )

                +     g(1)*g(1)/(2.*cw2) * (
          (-0.5 + sw2)     * ( cstau2* A0(mstau_1) + sstau2* A0(mstau_2) ) -
          sw2              * ( sstau2* A0(mstau_1) + cstau2* A0(mstau_2) ) +

          0.5 * A0(0) );

  /*
  // 1st 2nd family squarks and sleptons contribution
  for (int i=1; i<=2; i++)
    PI_22 = PI_22 + 3.* g(1)*g(1)/(2.*cw2) * (
            (0.5 - 2.*sw2/3.)*  A0(sqrt(real(mQu(i,i)))) +
            2.*sw2/3.        *  A0(sqrt(real(Mu2(i,i)))) +

            (-0.5 + sw2/3.)  *  A0(sqrt(real(mQd(i,i)))) -
            sw2/3.           *  A0(sqrt(real(Md2(i,i)))) )

                  +     g(1)*g(1)/(2.*cw2) * (
            (-0.5 + sw2)     *  A0(sqrt(real(ML2(i,i)))) -
            sw2              *  A0(sqrt(real(Me2(i,i)))) +

            0.5 * A0(mMv(i)) );
  */

  // 3rd line in Pierce ///////////////////////////////////////////////////////
  MatrixXd L_2LRu1(2,2), L_2LRu2(2,2), L_2LRu3(2,2),
                 L_2LRd1(2,2), L_2LRd2(2,2), L_2LRd3(2,2),
                 L_2LRe1(2,2), L_2LRe2(2,2), L_2LRe3(2,2);


  /*
  // Up squark : 1st generation
  L_2LRu1(1,1) = -g(1)*M_Z*sb*(0.5-2.*sw2/3.)/cw +sqrt(2.)*mYu(0)*mYu(0)/vu;
  L_2LRu1(1,2) = real(TYu(1,1))/sqrt(2.);
  L_2LRu1(2,1) = L_2LRu1(1,2);
  L_2LRu1(2,2) = -g(1)*M_Z*sb*2.*sw2/(3.*cw) + sqrt(2.)*mYu(0)*mYu(0)/vu;

    PI_22 = PI_22 + 3.* (
       L_2LRu1(1,1)*L_2LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
       L_2LRu1(2,2)*L_2LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))) );

  // Up squark : 2nd generation
  L_2LRu2(1,1) = -g(1)*M_Z*sb*(0.5-2.*sw2/3.)/cw +sqrt(2.)*mYu(1)*mYu(1)/vu;
  L_2LRu2(1,2) = real(TYu(2,2))/sqrt(2.);
  L_2LRu2(2,1) = L_2LRu2(1,2);
  L_2LRu2(2,2) = -g(1)*M_Z*sb*2.*sw2/(3.*cw) + sqrt(2.)*mYu(1)*mYu(1)/vu;

    PI_22 = PI_22 + 3.* (
       L_2LRu2(1,1)*L_2LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
       L_2LRu2(2,2)*L_2LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))) );
  */

  // Up squark : 3rd generation
  L_2LRu3(0,0) = -g(1)*M_Z*sb*(0.5-2.*sw2/3.)/cw +sqrt(2.)*mYu(2)*mYu(2)/vu;
  L_2LRu3(0,1) = real(TYu(2,2))/sqrt(2.);
  L_2LRu3(1,0) = L_2LRu3(0,1);
  L_2LRu3(1,1) = -g(1)*M_Z*sb*2.*sw2/(3.*cw) + sqrt(2.)*mYu(2)*mYu(2)/vu;

  // matrix in u1u2 basis
  L_2LRu3 = G_stop * L_2LRu3 * G_stop.transpose().eval();

  PI_22 = PI_22 + 3.* (
              L_2LRu3(0,0) * L_2LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_2LRu3(0,1) * L_2LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_2LRu3(1,1) * L_2LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st generation
  L_2LRd1(1,1) = -g(1)*M_Z*sb*(-0.5 + sw2/3.)/cw;
  L_2LRd1(1,2) = -mYd(0)*mu/(vd*sqrt(2.));
  L_2LRd1(2,1) = L_2LRd1(1,2);
  L_2LRd1(2,2) = g(1)*M_Z*sb*sw2/(3.*cw);

  PI_22 = PI_22 + 3.* (
       L_2LRd1(1,1)*L_2LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
       L_2LRd1(2,2)*L_2LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))) );

  // Down squark: 2nd generation
  L_2LRd2(1,1) = -g(1)*M_Z*sb*(-0.5 + sw2/3.)/cw;
  L_2LRd2(1,2) = -mYd(1)*mu/(vd*sqrt(2.));
  L_2LRd2(2,1) = L_2LRd2(1,2);
  L_2LRd2(2,2) = g(1)*M_Z*sb*sw2/(3.*cw);

  PI_22 = PI_22 + 3.* (
        L_2LRd2(1,1)*L_2LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
        L_2LRd2(2,2)*L_2LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))) );
  */

  // Down squark: 3rd generation
  L_2LRd3(0,0) = -g(1)*M_Z*sb*(-0.5 + sw2/3.)/cw;
  L_2LRd3(0,1) = -mYd(2)*mu/(vd*sqrt(2.));
  L_2LRd3(1,0) = L_2LRd3(0,1);
  L_2LRd3(1,1) = g(1)*M_Z*sb*sw2/(3.*cw);

  // matrix in u1u2 basis
  L_2LRd3 = (G_sbot * L_2LRd3 * G_sbot.transpose()).eval();

  PI_22 = PI_22 + 3.* (
              L_2LRd3(0,0) * L_2LRd3(0,0) * B0(p,msbot_1,msbot_1) +
          2.* L_2LRd3(0,1) * L_2LRd3(0,1) * B0(p,msbot_1,msbot_2) +
              L_2LRd3(1,1) * L_2LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Sleptons: 1st generation
  L_2LRe1(1,1) = -g(1)*M_Z*sb*(-0.5 + sw2)/cw;
  L_2LRe1(1,2) = -mYe(0)*mu/(vd*sqrt(2.));
  L_2LRe1(2,1) = L_2LRe1(1,2);
  L_2LRe1(2,2) = g(1)*M_Z*sb*sw2/cw;

  PI_22 = PI_22 + (
       L_2LRe1(1,1)*L_2LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
       L_2LRe1(2,2)*L_2LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))) );

  // Sleptons: 2nd generation
  L_2LRe2(1,1) = -g(1)*M_Z*sb*(-0.5 + sw2)/cw;
  L_2LRe2(1,2) = -mYe(1)*mu/(vd*sqrt(2.));
  L_2LRe2(2,1) = L_2LRe2(1,2);
  L_2LRe2(2,2) = g(1)*M_Z*sb*sw2/cw;

    PI_22 = PI_22 + (
       L_2LRe2(1,1)*L_2LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
       L_2LRe2(2,2)*L_2LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))) );
  */

  // Sleptons: 3rd generation
  L_2LRe3(0,0) = -g(1)*M_Z*sb*(-0.5 + sw2)/cw;
  L_2LRe3(0,1) = -mYe(2)*mu/(vd*sqrt(2.));
  L_2LRe3(1,0) = L_2LRe3(0,1);
  L_2LRe3(1,1) = g(1)*M_Z*sb*sw2/cw;

  // matrix in u1u2 basis
  L_2LRe3 = (G_stau * L_2LRe3 * G_stau.transpose()).eval();

  PI_22 = PI_22 + (
              L_2LRe3(0,0) * L_2LRe3(0,0) * B0(p,mstau_1,mstau_1) +
          2.* L_2LRe3(0,1) * L_2LRe3(0,1) * B0(p,mstau_1,mstau_2) +
              L_2LRe3(1,1) * L_2LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino
  for (int i=0; i<3; i++)
    PI_22 = PI_22 + (g(1)*M_Z*sb*(0.5)/cw)*(g(1)*M_Z*sb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_22 = PI_22 + 0.25 * g(1)*g(1) * (
          cb2* ( 2.* F(p,mHp,M_W) + F(p,mA,M_Z)/cw2 )+
          sb2* ( 2.* F(p,M_W,M_W) + F(p,M_Z,M_Z)/cw2 ) );

  PI_22 = PI_22 + (7./4.)*g(1)*g(1)*sb2* (
          2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 );

  PI_22 = PI_22 - g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2);

//*****************************************************************************8

  PI_22 = PI_22 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (sb*(3.*sa2-ca2) - cb*s2a)*(sb*(3.*sa2-ca2) - cb*s2a)*B0(p,mH,mH) +
          (sb*(3.*ca2-sa2) + cb*s2a)*(sb*(3.*ca2-sa2) + cb*s2a)*B0(p,mh,mh) +
          2.*(2.*sb*s2a - cb*c2a)*(2.*sb*s2a - cb*c2a)*B0(p,mH,mh)        +
          c2b*sb*c2b*sb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )                   +
          2.*s2b*sb*s2b*sb * B0(p,M_Z,mA) );

  PI_22 = PI_22 - 0.125*g(1)*g(1)/cw2 * (
          (3.*sa2-ca2) * A0(mH) +
          (3.*ca2-sa2) * A0(mh) -
          c2b * A0(M_Z) + c2b * A0(mA) );

//*****************************************************************************8

  PI_22 = PI_22 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          c2b*sb*c2b*sb * B0(p,M_W,M_W) +
          (c2b*sb+2.*cw2*sb)*(c2b*sb+2.*cw2*sb) * B0(p,mHp,mHp) +
          2.*(s2b*sb-cw2*cb)*(s2b*sb-cw2*cb) * B0(p,M_W,mHp) );

  PI_22 = PI_22 - 0.25*g(1)*g(1)/cw2 * (
          (cw2 - sw2*c2b) * A0(M_W) +
          (cw2 + sw2*c2b) * A0(mHp) );

//cout<<"PI_22 h= "<<PI_22<<endl;

//*****************************************************************************8
// neutralino contribution

  MatrixXcd Ann2(4,4), Bnn2(4,4);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      Ann2(i,j) = std::complex<double> (0.,0.);
      Bnn2(i,j) = std::complex<double> (0.,0.);
    }

  Ann2(0,3) = 0.5* g(0)*sqrt(3./5.);
  Ann2(3,0) = Ann2(0,3);
  Ann2(1,3) = -0.5* g(1);
  Ann2(3,1) = Ann2(1,3);

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      Bnn2(i,j) = Ann2(i,j);

// to mass eigenstates basis

  Ann2 = (Nt.conjugate() * Ann2 * Nt.adjoint()).eval();

  Bnn2 = (Nt * Bnn2 * Nt.transpose()).eval();

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_22 = PI_22 + 0.5*
              real( Ann2(i,j)*conj(Ann2(i,j)) + Bnn2(i,j)*conj(Bnn2(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann2(i,j)*conj(Bnn2(i,j)) +  conj(Ann2(i,j))*Bnn2(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//cout<<"PI_22 n= "<<PI_22<<endl;
//*****************************************************************************8
// chargino contribution

  MatrixXcd Acc2(2,2), Bcc2(2,2);

  Acc2(0,0) = std::complex<double> (0.,0.);
  Acc2(1,1) = std::complex<double> (0.,0.);
  Acc2(0,1) = std::complex<double> (0.,0.);
  Acc2(1,0) = g(1) / sqrt(2.) ;

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      Bcc2(i,j) =  Acc2(j,i);

// to mass eigenstates basis

  Acc2 = (V.conjugate() * Acc2 * U.adjoint()).eval();

  Bcc2 = U * Bcc2 * V.transpose().eval();

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_22 = PI_22 +
              std::real( Acc2(i,j)*conj(Acc2(i,j)) + Bcc2(i,j)*conj(Bcc2(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* std::real( Acc2(i,j)*conj(Bcc2(i,j)) + conj(Acc2(i,j))*Bcc2(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//cout<<"PI_22 c= "<<PI_22<<endl;
//*****************************************************************************8
// overall factor

  PI_22 = PI_22/(16.*pi*pi);

//  cout<<"PI_22 = "<<PI_22<<endl;


//*****************************************************************************8
// CP even Higgs self-energies 12

//  p = mh_pole;
//  p = M_Zpole;

  PI_12 = 0.;

  /*
  // Up squark : 1st 2nd generation
  PI_12 = PI_12 + 3.* (
     L_1LRu1(1,1)*L_2LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
     L_1LRu1(2,2)*L_2LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))) );

  PI_12 = PI_12 + 3.* (
     L_1LRu2(1,1)*L_2LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
     L_1LRu2(2,2)*L_2LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))) );
  */

  // Up squark : 3rd generation
  PI_12 = PI_12 + 3.* (
              L_1LRu3(0,0) * L_2LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_1LRu3(0,1) * L_2LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_1LRu3(1,1) * L_2LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st 2nd generation
  PI_12 = PI_12 + 3.* (
    L_1LRd1(1,1)*L_2LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
    L_1LRd1(2,2)*L_2LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))) );

  PI_12 = PI_12 + 3.* (
    L_1LRd2(1,1)*L_2LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
    L_1LRd2(2,2)*L_2LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))) );
  */

  // Down squark: 3rd generation
  PI_12 = PI_12 + 3.* (
                L_1LRd3(0,0) * L_2LRd3(0,0) * B0(p,msbot_1,msbot_1) +
            2.* L_1LRd3(0,1) * L_2LRd3(0,1) * B0(p,msbot_1,msbot_2) +
                L_1LRd3(1,1) * L_2LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Leptons: 1st 2nd generation
  PI_12 = PI_12 + (
    L_1LRe1(1,1)*L_2LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
    L_1LRe1(2,2)*L_2LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))) );

    PI_12 = PI_12 + (
    L_1LRe2(1,1)*L_2LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
    L_1LRe2(2,2)*L_2LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))) );
  */

  // Leptons: 3rd generation
  PI_12 = PI_12 + (
              L_1LRe3(0,0) * L_2LRe3(0,0) * B0(p,mstau_1,mstau_1) +
          2.* L_1LRe3(0,1) * L_2LRe3(0,1) * B0(p,mstau_1,mstau_2) +
              L_1LRe3(1,1) * L_2LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino
  for (int i=0; i<3; i++)
    PI_12 = PI_12 + (g(1)*M_Z*cb*(0.5)/cw)*(-g(1)*M_Z*sb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_12 = PI_12 + 0.25 * g(1)*g(1) * sb*cb * (
          2.* F(p,M_W,M_W) - 2.*F(p,mHp,M_W) +
          ( F(p,M_Z,M_Z) - F(p,mA,M_Z) )/cw2 +
          7.* ( 2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 ) );

//*****************************************************************************8

  PI_12 = PI_12 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (cb*(3.*ca2-sa2) - sb*s2a)*(sb*(3.*sa2-ca2) - cb*s2a)*B0(p,mH,mH) +
          (cb*(3.*sa2-ca2) + sb*s2a)*(sb*(3.*ca2-sa2) + cb*s2a)*B0(p,mh,mh) +
          2.*(-2.*cb*s2a - sb*c2a)*(2.*sb*s2a - cb*c2a)*B0(p,mH,mh)
          -c2b*cb*c2b*sb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )
          -2.*s2b*cb*s2b*sb * B0(p,M_Z,mA) );

  PI_12 = PI_12 - 0.125*g(1)*g(1)/cw2 *
           s2a * (A0(mh) - A0(mH));

//*****************************************************************************8

  PI_12 = PI_12 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          -c2b*cb*c2b*sb * B0(p,M_W,M_W) +
          (-c2b*cb+2.*cw2*cb)*(c2b*sb+2.*cw2*sb) * B0(p,mHp,mHp) +
          2.*(-s2b*cb+cw2*sb)*(s2b*sb-cw2*cb) * B0(p,M_W,mHp) );

  PI_12 = PI_12 - 0.25*g(1)*g(1)/cw2 *
          cw2 * s2b * ( - A0(M_W) + A0(mHp) );

//*****************************************************************************8
// neutralino contribution

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_12 = PI_12 + 0.5*
              real( Ann2(i,j)*conj(Ann1(i,j)) + Bnn2(i,j)*conj(Bnn1(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann2(i,j)*conj(Bnn1(i,j)) +  conj(Ann1(i,j))*Bnn2(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//*****************************************************************************8
// chargino contribution

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_12 = PI_12 +
              real( Acc2(i,j)*conj(Acc1(i,j)) + Bcc2(i,j)*conj(Bcc1(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* real( Acc2(i,j)*conj(Bcc1(i,j)) + conj(Acc1(i,j))*Bcc2(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//*****************************************************************************8
// overall factor

  PI_12 = PI_12/(16.*pi*pi);

//  cout<<"PI_12 = "<<PI_12<<endl;

// saving as PI_xxh

  PI_11h = PI_11;
  PI_22h = PI_22;
  PI_12h = PI_12;

//*****************************************************************************8
// CP even Higgs self-energies


// for heavy Higgs:

  p = mH_pole;

  PI_11 = 0.;

  // 1st line in Pierce ///////////////////////////////////////////////////////
  // Quarks and leptons contribution
  for (int i=0; i<3; i++)
    PI_11 = PI_11 +
            3.*mYd(i)*mYd(i)*( ( p*p-4.*mYd(i)*mYd(i) )* B0(p,mYd(i),mYd(i))
            + 2.* A0(mYd(i)) )/(vd*vd) +

            mYe(i)*mYe(i)*( ( p*p-4.*mYe(i)*mYe(i) )* B0(p,mYe(i),mYe(i))
            + 2.* A0(mYe(i)) )/(vd*vd);

  // 3rd family squarks and sleptons contribution
  PI_11 = PI_11 - 3.*mYd(2)*mYd(2)* ( A0(msbot_1) + A0(msbot_2) )/(vd*vd)
                -    mYe(2)*mYe(2)* ( A0(mstau_1) + A0(mstau_2) )/(vd*vd);

  // 2nd line in Pierce ///////////////////////////////////////////////////////
  // 3rd family squarks and sleptons contribution
  PI_11 = PI_11 - 3.* g(1)*g(1)/(2.*cw2) * (
          (0.5 - 2.*sw2/3.)* ( cstop2* A0(mstop_1) + sstop2* A0(mstop_2) ) +
          2.*sw2/3.        * ( sstop2* A0(mstop_1) + cstop2* A0(mstop_2) ) +

          (-0.5 + sw2/3.)  * ( csbot2* A0(msbot_1) + ssbot2* A0(msbot_2) ) -
          sw2/3.           * ( ssbot2* A0(msbot_1) + csbot2* A0(msbot_2) ) )

                -     g(1)*g(1)/(2.*cw2) * (
          (-0.5 + sw2)     * ( cstau2* A0(mstau_1) + sstau2* A0(mstau_2) ) -
          sw2              * ( sstau2* A0(mstau_1) + cstau2* A0(mstau_2) ) +

          0.5 * A0(0) );

  // 1st 2nd family squarks and sleptons contribution
  /*
  for (int i=1; i<=2; i++)
    PI_11 = PI_11 - 3.* g(1)*g(1)/(2.*cw2) * (
            (0.5 - 2.*sw2/3.)*  A0(sqrt(real(mQu(i,i)))) +
            2.*sw2/3.        *  A0(sqrt(real(Mu2(i,i)))) +

            (-0.5 + sw2/3.)  *  A0(sqrt(real(mQd(i,i)))) -
            sw2/3.           *  A0(sqrt(real(Md2(i,i)))) )

                  -     g(1)*g(1)/(2.*cw2) * (
            (-0.5 + sw2)     *  A0(sqrt(real(ML2(i,i)))) -
            sw2              *  A0(sqrt(real(Me2(i,i)))) +

            0.5 * A0(mMv(i)) );
  */

  // 3rd line in Pierce ///////////////////////////////////////////////////////
  /*
  // Up squark : 1st 2nd  generation
  PI_11 = PI_11 + 3.* (
     L_1LRu1(1,1)*L_1LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
     L_1LRu1(2,2)*L_1LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))) );

  PI_11 = PI_11 + 3.* (
     L_1LRu2(1,1)*L_1LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
     L_1LRu2(2,2)*L_1LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))) );
  */

  // Up squark : 3rd generation
  PI_11 = PI_11 + 3.* (
              L_1LRu3(0,0) * L_1LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_1LRu3(0,1) * L_1LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_1LRu3(1,1) * L_1LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st generation
  PI_11 = PI_11 + 3.* (
     L_1LRd1(1,1)*L_1LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
     L_1LRd1(2,2)*L_1LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))) );

  PI_11 = PI_11 + 3.* (
          L_1LRd2(1,1)*L_1LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
          L_1LRd2(2,2)*L_1LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))) );
  */

  // Down squark: 3rd generation
  PI_11 = PI_11 + 3.* (
              L_1LRd3(0,0) * L_1LRd3(0,0) * B0(p,msbot_1,msbot_1) +
          2.* L_1LRd3(0,1) * L_1LRd3(0,1) * B0(p,msbot_1,msbot_2) +
              L_1LRd3(1,1) * L_1LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Sleptons: 1st generation
  PI_11 = PI_11 + (
     L_1LRe1(1,1)*L_1LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
     L_1LRe1(2,2)*L_1LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))) );

  PI_11 = PI_11 + (
     L_1LRe2(1,1)*L_1LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
     L_1LRe2(2,2)*L_1LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))) );
  */

  // Sleptons: 3rd generation
  PI_11 = PI_11 + (
              L_1LRe3(0,0) * L_1LRe3(0,0) * B0(p,mstau_1,mstau_1) +
          2.* L_1LRe3(0,1) * L_1LRe3(0,1) * B0(p,mstau_1,mstau_2) +
              L_1LRe3(1,1) * L_1LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino contribution
  for (int i=0; i<3; i++)
    PI_11 = PI_11 + (g(1)*M_Z*cb*(0.5)/cw)*(g(1)*M_Z*cb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_11 = PI_11 + 0.25 * g(1)*g(1) * (
          sb2* ( 2.* F(p,mHp,M_W) + F(p,mA,M_Z)/cw2 )+
          cb2* ( 2.* F(p,M_W,M_W) + F(p,M_Z,M_Z)/cw2 ) );

  PI_11 = PI_11 + (7./4.)*g(1)*g(1)*cb2* (
          2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 );

  PI_11 = PI_11 - g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2);

//*****************************************************************************8

  PI_11 = PI_11 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (cb*(3.*ca2-sa2) - sb*s2a)*(cb*(3.*ca2-sa2) - sb*s2a)*B0(p,mH,mH) +
          (cb*(3.*sa2-ca2) + sb*s2a)*(cb*(3.*sa2-ca2) + sb*s2a)*B0(p,mh,mh) +
          2.*(-2.*cb*s2a - sb*c2a)*(-2.*cb*s2a - sb*c2a)*B0(p,mH,mh)        +
          c2b*cb*c2b*cb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )                   +
          2.*s2b*cb*s2b*cb * B0(p,M_Z,mA) );

  PI_11 = PI_11 - 0.125*g(1)*g(1)/cw2 * (
          (3.*ca2-sa2) * A0(mH) +
          (3.*sa2-ca2) * A0(mh) +
          c2b * A0(M_Z) - c2b * A0(mA) );

//*****************************************************************************8

  PI_11 = PI_11 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          c2b*cb*c2b*cb * B0(p,M_W,M_W) +
          (-c2b*cb+2.*cw2*cb)*(-c2b*cb+2.*cw2*cb) * B0(p,mHp,mHp) +
          2.*(-s2b*cb+cw2*sb)*(-s2b*cb+cw2*sb) * B0(p,M_W,mHp) );

  PI_11 = PI_11 - 0.25*g(1)*g(1)/cw2 * (
          (cw2 + sw2*c2b) * A0(M_W) +
          (cw2 - sw2*c2b) * A0(mHp) );

//*****************************************************************************8
// neutralino contribution

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_11 = PI_11 + 0.5*
              real( Ann1(i,j)*conj(Ann1(i,j)) + Bnn1(i,j)*conj(Bnn1(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann1(i,j)*conj(Bnn1(i,j)) +  conj(Ann1(i,j))*Bnn1(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//*****************************************************************************8
// chargino contribution

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_11 = PI_11 +
              real( Acc1(i,j)*conj(Acc1(i,j)) + Bcc1(i,j)*conj(Bcc1(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* real( Acc1(i,j)*conj(Bcc1(i,j)) + conj(Acc1(i,j))*Bcc1(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//*****************************************************************************8
// overall factor

   PI_11 = PI_11/(16.*pi*pi);

//cout<<"PI_11 = "<<PI_11<<endl;

//*****************************************************************************8
// CP even Higgs self-energies 22

  PI_22 = 0.;

  // 1st line in Pierce ///////////////////////////////////////////////////////
  // Quarks and leptons contribution
  for (int i=0; i<3; i++)
    PI_22 = PI_22 +
            3.*mYu(i)*mYu(i)*( ( p*p-4.*mYu(i)*mYu(i) )* B0(p,mYu(i),mYu(i))
            + 2.* A0(mYu(i)) )/(vu*vu) ;

  // 3rd family squarks and sleptons contribution
  PI_22 = PI_22 - 3.*mYu(2)*mYu(2)* ( A0(mstop_1) + A0(mstop_2) )/(vu*vu);

  // 2nd line in Pierce ///////////////////////////////////////////////////////
  // 3rd family squarks and sleptons contribution0
  PI_22 = PI_22 + 3.* g(1)*g(1)/(2.*cw2) * (
          (0.5 - 2.*sw2/3.)* ( cstop2* A0(mstop_1) + sstop2* A0(mstop_2) ) +
          2.*sw2/3.        * ( sstop2* A0(mstop_1) + cstop2* A0(mstop_2) ) +

          (-0.5 + sw2/3.)  * ( csbot2* A0(msbot_1) + ssbot2* A0(msbot_2) ) -
          sw2/3.           * ( ssbot2* A0(msbot_1) + csbot2* A0(msbot_2) ) )

                +     g(1)*g(1)/(2.*cw2) * (
          (-0.5 + sw2)     * ( cstau2* A0(mstau_1) + sstau2* A0(mstau_2) ) -
          sw2              * ( sstau2* A0(mstau_1) + cstau2* A0(mstau_2) ) +

          0.5 * A0(0) );

  /*
  // 1st 2nd family squarks and sleptons contribution
  for (int i=1; i<=2; i++)
    PI_22 = PI_22 + 3.* g(1)*g(1)/(2.*cw2) * (
            (0.5 - 2.*sw2/3.)*  A0(sqrt(real(mQu(i,i)))) +
            2.*sw2/3.        *  A0(sqrt(real(Mu2(i,i)))) +

            (-0.5 + sw2/3.)  *  A0(sqrt(real(mQd(i,i)))) -
            sw2/3.           *  A0(sqrt(real(Md2(i,i)))) )

            + g(1)*g(1)/(2.*cw2) * (
            (-0.5 + sw2)     *  A0(sqrt(real(ML2(i,i)))) -
            sw2              *  A0(sqrt(real(Me2(i,i)))) +

            0.5 * A0(mMv(i)) );
  */

  // 3rd line in Pierce ///////////////////////////////////////////////////////
  /*
  // Up squark : 1st 2nd generation
  PI_22 = PI_22 + 3.* (
    L_2LRu1(1,1)*L_2LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
    L_2LRu1(2,2)*L_2LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))));

  PI_22 = PI_22 + 3.* (
    L_2LRu2(1,1)*L_2LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
    L_2LRu2(2,2)*L_2LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))));
  */

  // Up squark : 3rd generation
  PI_22 = PI_22 + 3.* (
              L_2LRu3(0,0) * L_2LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_2LRu3(0,1) * L_2LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_2LRu3(1,1) * L_2LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st 2nd generation
  PI_22 = PI_22 + 3.* (
    L_2LRd1(1,1)*L_2LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
    L_2LRd1(2,2)*L_2LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))));
  PI_22 = PI_22 + 3.* (
    L_2LRd2(1,1)*L_2LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
    L_2LRd2(2,2)*L_2LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))));
  */

  // Down squark: 3rd generation
  PI_22 = PI_22 + 3.* (
              L_2LRd3(0,0) * L_2LRd3(0,0) * B0(p,msbot_1,msbot_1) +
          2.* L_2LRd3(0,1) * L_2LRd3(0,1) * B0(p,msbot_1,msbot_2) +
              L_2LRd3(1,1) * L_2LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Sleptons: 1st 2nd generation
  PI_22 = PI_22 + (
    L_2LRe1(1,1)*L_2LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
    L_2LRe1(2,2)*L_2LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))));
  PI_22 = PI_22 + (
    L_2LRe2(1,1)*L_2LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
    L_2LRe2(2,2)*L_2LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))));
  */

  // Sleptons: 3rd generation
  PI_22 = PI_22 + (
              L_2LRe3(0,0) * L_2LRe3(0,0) * B0(p,mstau_1,mstau_1) +
          2.* L_2LRe3(0,1) * L_2LRe3(0,1) * B0(p,mstau_1,mstau_2) +
              L_2LRe3(1,1) * L_2LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino
  for (int i=0; i<3; i++)
    PI_22 = PI_22 + (g(1)*M_Z*sb*(0.5)/cw)*(g(1)*M_Z*sb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_22 = PI_22 + 0.25 * g(1)*g(1) * (
          cb2* ( 2.* F(p,mHp,M_W) + F(p,mA,M_Z)/cw2 )+
          sb2* ( 2.* F(p,M_W,M_W) + F(p,M_Z,M_Z)/cw2 ) );

  PI_22 = PI_22 + (7./4.)*g(1)*g(1)*sb2* (
          2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 );

  PI_22 = PI_22 - g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2);

  PI_22 = PI_22 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (sb*(3.*sa2-ca2) - cb*s2a)*(sb*(3.*sa2-ca2) - cb*s2a)*B0(p,mH,mH) +
          (sb*(3.*ca2-sa2) + cb*s2a)*(sb*(3.*ca2-sa2) + cb*s2a)*B0(p,mh,mh) +
          2.*(2.*sb*s2a - cb*c2a)*(2.*sb*s2a - cb*c2a)*B0(p,mH,mh)        +
          c2b*sb*c2b*sb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )                   +
          2.*s2b*sb*s2b*sb * B0(p,M_Z,mA) );

  PI_22 = PI_22 - 0.125*g(1)*g(1)/cw2 * (
          (3.*sa2-ca2) * A0(mH) +
          (3.*ca2-sa2) * A0(mh) -
          c2b * A0(M_Z) + c2b * A0(mA) );

  PI_22 = PI_22 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          c2b*sb*c2b*sb * B0(p,M_W,M_W) +
          (c2b*sb+2.*cw2*sb)*(c2b*sb+2.*cw2*sb) * B0(p,mHp,mHp) +
          2.*(s2b*sb-cw2*cb)*(s2b*sb-cw2*cb) * B0(p,M_W,mHp) );

  PI_22 = PI_22 - 0.25*g(1)*g(1)/cw2 * (
          (cw2 - sw2*c2b) * A0(M_W) +
          (cw2 + sw2*c2b) * A0(mHp) );

//*****************************************************************************8
// neutralino contribution

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_22 = PI_22 + 0.5*
              real( Ann2(i,j)*conj(Ann2(i,j)) + Bnn2(i,j)*conj(Bnn2(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann2(i,j)*conj(Bnn2(i,j)) +  conj(Ann2(i,j))*Bnn2(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//*****************************************************************************8
// chargino contribution

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_22 = PI_22 +
              real( Acc2(i,j)*conj(Acc2(i,j)) + Bcc2(i,j)*conj(Bcc2(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* real( Acc2(i,j)*conj(Bcc2(i,j)) + conj(Acc2(i,j))*Bcc2(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

//*****************************************************************************8
// overall factor

  PI_22 = PI_22/(16.*pi*pi);

//  cout<<"PI_22 = "<<PI_22<<endl;


//*****************************************************************************8
// CP even Higgs self-energies 12

  PI_12 = 0.;

  /*
  // Up squark : 1st 2nd generation
  PI_12 = PI_12 + 3.* (
     L_1LRu1(1,1)*L_2LRu1(1,1) * B0(p,sqrt(real(mQu(1,1))),sqrt(real(mQu(1,1)))) +
     L_1LRu1(2,2)*L_2LRu1(2,2) * B0(p,sqrt(real(Mu2(1,1))),sqrt(real(Mu2(1,1)))) );

  PI_12 = PI_12 + 3.* (
     L_1LRu2(1,1)*L_2LRu2(1,1) * B0(p,sqrt(real(mQu(2,2))),sqrt(real(mQu(2,2)))) +
     L_1LRu2(2,2)*L_2LRu2(2,2) * B0(p,sqrt(real(Mu2(2,2))),sqrt(real(Mu2(2,2)))) );
  */

  // Up squark : 3rd generation
  PI_12 = PI_12 + 3.* (
              L_1LRu3(0,0) * L_2LRu3(0,0) * B0(p,mstop_1,mstop_1) +
          2.* L_1LRu3(0,1) * L_2LRu3(0,1) * B0(p,mstop_1,mstop_2) +
              L_1LRu3(1,1) * L_2LRu3(1,1) * B0(p,mstop_2,mstop_2) );

  /*
  // Down squark: 1st 2nd generation
  PI_12 = PI_12 + 3.* (
       L_1LRd1(1,1)*L_2LRd1(1,1) * B0(p,sqrt(real(mQd(1,1))),sqrt(real(mQd(1,1)))) +
       L_1LRd1(2,2)*L_2LRd1(2,2) * B0(p,sqrt(real(Md2(1,1))),sqrt(real(Md2(1,1)))) );

  PI_12 = PI_12 + 3.* (
        L_1LRd2(1,1)*L_2LRd2(1,1) * B0(p,sqrt(real(mQd(2,2))),sqrt(real(mQd(2,2)))) +
        L_1LRd2(2,2)*L_2LRd2(2,2) * B0(p,sqrt(real(Md2(2,2))),sqrt(real(Md2(2,2)))) );
  */

  // Down squark: 3rd generation
  PI_12 = PI_12 + 3.* (
                L_1LRd3(0,0) * L_2LRd3(0,0) * B0(p,msbot_1,msbot_1) +
            2.* L_1LRd3(0,1) * L_2LRd3(0,1) * B0(p,msbot_1,msbot_2) +
               L_1LRd3(1,1) * L_2LRd3(1,1) * B0(p,msbot_2,msbot_2) );

  /*
  // Leptons: 1st 2nd generation
  PI_12 = PI_12 + (
       L_1LRe1(1,1)*L_2LRe1(1,1) * B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) +
       L_1LRe1(2,2)*L_2LRe1(2,2) * B0(p,sqrt(real(Me2(1,1))),sqrt(real(Me2(1,1)))) );

    PI_12 = PI_12 + (
       L_1LRe2(1,1)*L_2LRe2(1,1) * B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) +
       L_1LRe2(2,2)*L_2LRe2(2,2) * B0(p,sqrt(real(Me2(2,2))),sqrt(real(Me2(2,2)))) );
  */

  // Leptons: 3rd generation
    PI_12 = PI_12 + (
                L_1LRe3(0,0) * L_2LRe3(0,0) * B0(p,mstau_1,mstau_1) +
            2.* L_1LRe3(0,1) * L_2LRe3(0,1) * B0(p,mstau_1,mstau_2) +
                L_1LRe3(1,1) * L_2LRe3(1,1) * B0(p,mstau_2,mstau_2) );

  // Sneutrino
  for (int i=0; i<3; i++)
    PI_12 = PI_12 + (g(1)*M_Z*cb*(0.5)/cw)*(-g(1)*M_Z*sb*(0.5)/cw)*
                     B0(p,0,0) ;

//*****************************************************************************8

  PI_12 = PI_12 + 0.25 * g(1)*g(1) * sb*cb * (
          2.* F(p,M_W,M_W) - 2.*F(p,mHp,M_W) +
          ( F(p,M_Z,M_Z) - F(p,mA,M_Z) )/cw2 +
          7.* ( 2.*M_W2* B0(p,M_W,M_W) + M_Z2* B0(p,M_Z,M_Z)/cw2 ) );

//*****************************************************************************8

  PI_12 = PI_12 + 0.125*g(1)*g(1)*M_Z2/cw2 * (
          (cb*(3.*ca2-sa2) - sb*s2a)*(sb*(3.*sa2-ca2) - cb*s2a)*B0(p,mH,mH) +
          (cb*(3.*sa2-ca2) + sb*s2a)*(sb*(3.*ca2-sa2) + cb*s2a)*B0(p,mh,mh) +
          2.*(-2.*cb*s2a - sb*c2a)*(2.*sb*s2a - cb*c2a)*B0(p,mH,mh)
          -c2b*cb*c2b*sb * ( B0(p,M_Z,M_Z) + B0(p,mA,mA) )
          -2.*s2b*cb*s2b*sb * B0(p,M_Z,mA) );

  PI_12 = PI_12 - 0.125*g(1)*g(1)/cw2 *
           s2a * (A0(mh) - A0(mH));

//*****************************************************************************8

  PI_12 = PI_12 + 0.25*g(1)*g(1)*M_Z2/cw2 * (
          -c2b*cb*c2b*sb * B0(p,M_W,M_W) +
         (-c2b*cb+2.*cw2*cb)*(c2b*sb+2.*cw2*sb) * B0(p,mHp,mHp) +
          2.*(-s2b*cb+cw2*sb)*(s2b*sb-cw2*cb) * B0(p,M_W,mHp) );

  PI_12 = PI_12 - 0.25*g(1)*g(1)/cw2 *
          cw2 * s2b * ( - A0(M_W) + A0(mHp) );


//*****************************************************************************8
// neutralino contribution

  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
      PI_12 = PI_12 + 0.5*
              real( Ann2(i,j)*conj(Ann1(i,j)) + Bnn2(i,j)*conj(Bnn1(i,j)) )
              * G(p,mMnt(i),mMnt(j)) -
              real( Ann2(i,j)*conj(Bnn1(i,j)) +  conj(Ann1(i,j))*Bnn2(i,j) )
              * mMnt(i)*mMnt(j) * B0(p,mMnt(i),mMnt(j));

//*****************************************************************************8
// chargino contribution

  for (int i=0; i<2; i++)
    for (int j=0; j<2; j++)
      PI_12 = PI_12 +
              real( Acc2(i,j)*conj(Acc1(i,j)) + Bcc2(i,j)*conj(Bcc1(i,j)) )
              * G(p,mMch(i),mMch(j)) -
              2.* real( Acc2(i,j)*conj(Bcc1(i,j)) + conj(Acc1(i,j))*Bcc2(i,j) )
              * mMch(i)*mMch(j) * B0(p,mMch(i),mMch(j));

  // overall factor
  PI_12 = PI_12/(16.*pi*pi);

  /////////////////////////////////////////////////////////////////////////////
  // tad poles t1/v1, t2/v2 : Pierce et al ////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // opposite overall sign due to - in A0

  // t1_v1 ////////////////////////////////////////////////////////////////////
  t1_v1 = 0.;

  for (int i=0; i<3; i++) {
    t1_v1 = t1_v1 + 6.* A0(mYd(i)) * mYd(i)*mYd(i)/(vd*vd)
                  + 2.* A0(mYe(i)) * mYe(i)*mYe(i)/(vd*vd);
  }

  // squarks and sleptons : 3rd generation mixing only approximation
   t1_v1 = t1_v1 - 3.*g(1)/(2.*M_W*cb) * A0(mstop_1) *
              ( cstop2*g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw
              - s2stop * mu*mYu(2)/(std::sqrt(2.)*vu)
              + sstop2*g(1)*M_Z*cb* 2.*sw2/(3.*cw) )

              - 3.*g(1)/(2.*M_W*cb) * A0(mstop_2) *
              ( sstop2*g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw
              + s2stop * mu*mYu(2)/(std::sqrt(2.)*vu)
              + cstop2*g(1)*M_Z*cb* 2.*sw2/(3.*cw) )

              - 3.*g(1)/(2.*M_W*cb) * A0(msbot_1) *
              ( csbot2* ( g(1)*M_Z*cb*(-0.5 + sw2/3.)/cw
                        + std::sqrt(2.)*mYd(2)*mYd(2)/vd )
              + s2sbot * real(TYd(2,2))/std::sqrt(2.)
              + ssbot2* ( -g(1)*M_Z*cb* sw2/(3.*cw)
                        + std::sqrt(2.)*mYd(2)*mYd(2)/vd ) )

              - 3.*g(1)/(2.*M_W*cb) * A0(msbot_2) *
              ( ssbot2* ( g(1)*M_Z*cb*(-0.5 + sw2/3.)/cw
                        + std::sqrt(2.)*mYd(2)*mYd(2)/vd )
              - s2sbot * real(TYd(2,2))/std::sqrt(2.)
              + csbot2* ( - g(1)*M_Z*cb* sw2/(3.*cw)
                        + std::sqrt(2.)*mYd(2)*mYd(2)/vd ) )

              - g(1)/(2.*M_W*cb) * A0(mstau_1) *
              ( cstau2* ( g(1)*M_Z*cb*(-0.5 + sw2)/cw
                        + std::sqrt(2.)*mYe(2)*mYe(2)/vd )
              + s2stau * real(TYe(2,2))/std::sqrt(2.)
              + sstau2* ( - g(1)*M_Z*cb* sw2/cw
                       + std::sqrt(2.)*mYe(2)*mYe(2)/vd ) )

              - g(1)/(2.*M_W*cb) * A0(mstau_2) *
              ( sstau2* ( g(1)*M_Z*cb*(-0.5 + sw2)/cw
                        + std::sqrt(2.)*mYe(2)*mYe(2)/vd )
              - s2stau * real(TYe(2,2))/std::sqrt(2.)
              + cstau2* ( - g(1)*M_Z*cb* sw2/cw
                        + std::sqrt(2.)*mYe(2)*mYe(2)/vd ) );

  // squarks and sleptons : 1st 2nd generation no mixing approximation
  // Need to keep this because this is the correction to the Higgs potential
   for (int i=0; i<2; i++) {
      t1_v1 = t1_v1 - 3.*g(1)/(2.*M_W*cb) *
                        ( A0(std::sqrt(std::real(mQu(i,i)))) *
                            g(1)*M_Z*cb*(0.5 - 2.*sw2/3.)/cw
                        + A0(std::sqrt(std::real(Mu2(i,i)))) *
                            g(1)*M_Z*cb* 2.*sw2/(3.*cw)
                        )
                    - 3.*g(1)/(2.*M_W*cb) *
                      ( A0(std::sqrt(std::real(mQd(i,i)))) *
                          ( g(1)*M_Z*cb*(-0.5 + sw2/3.)/cw
                          + std::sqrt(2.)*mYd(i)*mYd(i)/vd )
                      - A0(std::sqrt(std::real(Md2(i,i)))) *
                          ( g(1)*M_Z*cb* sw2/(3.*cw)
                          + std::sqrt(2.)*mYd(i)*mYd(i)/vd )
                      )
                    - g(1)/(2.*M_W*cb) *
                      ( A0(std::sqrt(std::real(ML2(i,i)))) *
                          ( g(1)*M_Z*cb*(-0.5 + sw2)/cw
                          + std::sqrt(2.)*mYe(i)*mYe(i)/vd )
                      - A0(std::sqrt(std::real(Me2(i,i)))) *
                          ( g(1)*M_Z*cb* sw2/cw
                          + std::sqrt(2.)*mYe(i)*mYe(i)/vd )
                      );
  }

  // Charged Higgs contribution
  t1_v1 = t1_v1 + g(1)*g(1)*c2b * (A0(mA) + 2.*A0(mHp)) / (8.*cw2)
          - g(1)*g(1) * A0(mHp)/2.
          - g(1)*g(1)* (3.*sa2 - ca2 + s2a*tanb) * A0(mh) / (8.*cw2)
          - g(1)*g(1)* (3.*ca2 - sa2 - s2a*tanb) * A0(mH) / (8.*cw2);

  // Neutrino contribution
  for (int i=0; i<4; i++) {
    t1_v1 = t1_v1 + g(1)*g(1) * mMnt(i) * A0(mMnt(i)) * real( Nt(i,2)*(Nt(i,1) - Nt(i,0)*sw/cw) )/(M_W*cb);
  }

  // Chargino contribution
  for (int i=0; i<2; i++) {
    t1_v1 = t1_v1 + sqrt(2.)*g(1)*g(1) * mMch(i) * A0(mMch(i)) * real( V(i,0)*U(i,1) )/(M_W*cb);
  }

  // W & Z contribution
  t1_v1 = t1_v1 - 0.75* g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2)
          - g(1)*g(1)*c2b * (2.*A0(M_W) + A0(M_Z)) /(8.*cw2);

  t1_v1 = t1_v1/(16.*pi*pi);


  // t2_v2 ////////////////////////////////////////////////////////////////////
  t2_v2 = 0.;

  for (int i=0; i<3; i++)
    t2_v2 = t2_v2 + 6.* A0(mYu(i)) * mYu(i)*mYu(i)/(vu*vu);

  // squarks and sleptons : 3rd generation mixing only approximation
  t2_v2 = t2_v2 - 3.*g(1)/(2.*M_W*sb) * A0(mstop_1) *
              ( cstop2*(-g(1)*M_Z*sb*(0.5 - 2.*sw2/3.)/cw
              + sqrt(2.)*mYu(2)*mYu(2)/vu )
              + s2stop* real(TYu(2,2))/sqrt(2.)
              + sstop2* ( - g(1)*M_Z*sb* 2.*sw2/(3.*cw)
              + sqrt(2.)*mYu(2)*mYu(2)/vu ) )

              - 3.*g(1)/(2.*M_W*sb) * A0(mstop_2) *
              ( sstop2*(-g(1)*M_Z*sb*(0.5 - 2.*sw2/3.)/cw
              + sqrt(2.)*mYu(2)*mYu(2)/vu )
              - s2stop* real(TYu(2,2))/sqrt(2.)
              + cstop2* ( - g(1)*M_Z*sb* 2.*sw2/(3.*cw)
              + sqrt(2.)*mYu(2)*mYu(2)/vu ) )


              - 3.*g(1)/(2.*M_W*sb) * A0(msbot_1) *
              ( csbot2*(-1.)* g(1)*M_Z*sb*( -0.5 + sw2/3.)/cw
              - s2sbot * mYd(2)*mu/(vd*sqrt(2.))
              + ssbot2*g(1)*M_Z*sb* sw2/(3.*cw) )

             - 3.*g(1)/(2.*M_W*sb) * A0(msbot_2) *
              ( ssbot2*(-1.)* g(1)*M_Z*sb*( -0.5 + sw2/3.)/cw
              + s2sbot * mYd(2)*mu/(vd*sqrt(2.))
              + csbot2*g(1)*M_Z*sb* sw2/(3.*cw) )


              - g(1)/(2.*M_W*sb) * A0(mstau_1) *
              ( cstau2*(-1.)* g(1)*M_Z*sb*(-0.5 + sw2)/cw
              -  s2stau *mYe(2)*mu/(vd*sqrt(2.))
              + sstau2*g(1)*M_Z*sb* sw2/cw )

              - g(1)/(2.*M_W*sb) * A0(mstau_2) *
              ( sstau2*(-1.)* g(1)*M_Z*sb*(-0.5 + sw2)/cw
              +  s2stau *mYe(2)*mu/(vd*sqrt(2.))
              + cstau2*g(1)*M_Z*sb* sw2/cw );

  // squarks and sleptons : 1st 2nd generation no mixing approximation
  // Need to keep this because this is the correction to the Higgs potential
  for (int i=0; i<2; i++) {
    t2_v2 = t2_v2 -  3.*g(1)/(2.*M_W*sb) *
                      ( A0(sqrt(real(mQu(i,i)))) *
                          (-g(1)*M_Z*sb*(0.5-2.*sw2/3.)/cw
                          + sqrt(2.)*mYu(i)*mYu(i)/vu)
                      + A0(sqrt(real(Mu2(i,i)))) *
                          (-g(1)*M_Z*sb*2.*sw2/(3.*cw)
                          + sqrt(2.)*mYu(i)*mYu(i)/vu)
                      )
                  - 3.*g(1)/(2.*M_W*sb)  *
                      ( A0(sqrt(real(mQd(i,i)))) *
                        (-1.)* g(1)*M_Z*sb*( -0.5 + sw2/3.)/cw
                      + A0(sqrt(real(Md2(i,i)))) *
                        g(1)*M_Z*sb* sw2/(3.*cw) )
                  - g(1)/(2.*M_W*sb) *
                      ( A0(sqrt(real(ML2(i,i)))) *
                        (-1.)* g(1)*M_Z*sb*(-0.5 + sw2)/cw
                      + A0(sqrt(real(Me2(i,i)))) *
                       g(1)*M_Z*sb* sw2/cw );
  }

  for (int a=0; a<3; a++)
    t2_v2 = t2_v2 - g(1)/(2.*M_W*sb) * A0(0)
            * (-1.)*g(1)*M_Z*sb*0.5/cw;

  t2_v2 = t2_v2 - g(1)*g(1)*c2b * (A0(mA) + 2.*A0(mHp)) / (8.*cw2)
          - g(1)*g(1) * A0(mHp)/2.
          - g(1)*g(1)* (3.*ca2 - sa2 + s2a/tanb) * A0(mh) / (8.*cw2)
          - g(1)*g(1)* (3.*sa2 - ca2 - s2a/tanb) * A0(mH) / (8.*cw2);

  for (int i=0; i<4; i++)
    t2_v2 = t2_v2 - g(1)*g(1) * mMnt(i) * A0(mMnt(i)) *
            real( Nt(i,3)*(Nt(i,1) - Nt(i,0)*sw/cw) )/(M_W*sb);

  for (int i=0; i<2; i++)
    t2_v2 = t2_v2 + sqrt(2.)*g(1)*g(1) * mMch(i) * A0(mMch(i)) *
            real( V(i,1)*U(i,0) )/(M_W*sb);

  t2_v2 = t2_v2 - 0.75* g(1)*g(1)* (2.*A0(M_W) + A0(M_Z)/cw2)
          + g(1)*g(1)*c2b * (2.*A0(M_W) + A0(M_Z)) /(8.*cw2);

  t2_v2 = t2_v2/(16.*pi*pi);

  /////////////////////////////////////////////////////////////////////////////
  // EWSB Consistency Check ///////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  // mHu2bar, mHd2bar
  mHd2bar = mHd - t1_v1;
  mHu2bar = mHu - t2_v2;

  // EWSB condition and equation for m_A
  mju2_new = 0.5*( t2b*(mHu2bar*tanb - mHd2bar/tanb) - M_Z2 );
  m_A2_new = (mHu2bar - mHd2bar)/c2b - M_Z2;

  double mA2_pole;
  bA = sb2*t1_v1 + cb2*t2_v2;
  mA2_pole = m_A2_new  - PI_AA + bA;
  if (mA2_pole >= 0.) {
    mA_pole = sqrt(mA2_pole);
  } else {
    std::cout<<"WARNING!!! Negative mA_pole^2"<<std::endl;
    Pe = Pe*10.;
    mA_pole = sqrt(-mA2_pole);
    return Pe;
  }

  mHd_new  = sb2*(mA2 + PI_AA + M_Z2 - bA) - 0.5*(2.*mu*mu + M_Z2) + t1_v1;
  mHu_new  = cb2*(mA2 + PI_AA + M_Z2 - bA) - 0.5*(2.*mu*mu + M_Z2) + t2_v2;

  /////////////////////////////////////////////////////////////////////////////
  // new Higgses and mixing ///////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  if (m_A2_new >0.)
  {
    //Tree level Higgses mass (Calc from corrected tree-level mA)
    mHp2  = m_A2_new + M_W2;
    mHp  = sqrt(mHp2);

    mH2  = 0.5*(m_A2_new + M_Z2 +
           sqrt((m_A2_new+M_Z2)*(m_A2_new+M_Z2)
           - 4.*m_A2_new*M_Z2*c2b*c2b));
    mH   = sqrt(mH2);

    mh2  = 0.5*(m_A2_new + M_Z2 -
           sqrt((m_A2_new+M_Z2)*(m_A2_new+M_Z2)
           - 4.*m_A2_new*M_Z2*c2b*c2b));
    mh   = sqrt(mh2);

    //1-loop Higgses mass
    t2a   = t2b * (m_A2_new + M_Z2)/(m_A2_new - M_Z2);
    c2a   = - c2b*(m_A2_new - M_Z2)/(mH2 - mh2);
    s2a   = t2a * c2a;

    double a, b, d;
    double mH_pole2, mh_pole2, mHp_pole2;

    //Heavy Higgs
    a = M_Z2*cb2 + mA2*sb2 - PI_11 + t1_v1;
    b = - (M_Z2 + mA2)*sb*cb - PI_12;
    d = M_Z2*sb2 + mA2*cb2 - PI_22 + t2_v2;

    mH_pole2  = 0.5*( a+d + sqrt((a-d)*(a-d) + 4.*b*b) );
    if ( mH_pole2 > 0.) mH_pole   = sqrt(mH_pole2);
    else mH_pole   = 0.;

    //Light Higgs
    a = M_Z2*cb2 + mA2*sb2 - PI_11h + t1_v1;
    b = - (M_Z2 + mA2)*sb*cb - PI_12h;
    d = M_Z2*sb2 + mA2*cb2 - PI_22h + t2_v2;

    mh_pole2  = 0.5*( a+d - sqrt((a-d)*(a-d) + 4.*b*b) );
    if ( mh_pole2 > 0.) mh_pole   = sqrt(mh_pole2);
    else mh_pole   = 0.;

    //Charged Higgs (Approx: Calc does not include PI_HH, PI_AA, PI_WW)
    mHp_pole2 = mA2_pole + M_Wpole2;
    if (mHp_pole2 > 0.) mHp_pole = sqrt(mHp_pole2);
    else mHp_pole = 0.;
  }


  ///////////////////////////////////////////////////////////////////////////
  // Integrating out third family scalars and gluinos ///////////////////////
  ///////////////////////////////////////////////////////////////////////////
  double invalpha1_un, invalpha1_c, invalpha2_un, invalpha2_c, invalpha3_un, invalpha3_c;

  invalpha1_un = (4.*pi)/(g(0)*g(0));
  invalpha2_un = (4.*pi)/(g(1)*g(1));
  invalpha3_un = (4.*pi)/(g(2)*g(2));
  // Ignored that for stop mass and gauge eigenstates are not the same
  invalpha3_c = invalpha3_un - 1./(12.*pi)*log(mMu(4)*mMu(5)*mMd(4)*mMd(5)/pow(MEW,4))- 1./pi*log(Mg/MEW);
  invalpha2_c = invalpha2_un - 1./(4.*pi)*log(mMu(4)/MEW)- 1./(12.*pi)*log(mMe(4)/MEW);
  invalpha1_c = invalpha1_un - 1./(60.*pi)*log(mMu(4)/MEW)-  2./(15.*pi)*log(mMu(5)/MEW)- \
    1./(30.*pi)*log(mMd(5)/MEW)- 1./(20.*pi)*log(mMe(4)/MEW)- 1./(10.*pi)*log(mMe(5)/MEW);

  g(0) = sqrt(4.*pi/invalpha1_c);
  g(1) = sqrt(4.*pi/invalpha2_c);
  g(2) = sqrt(4.*pi/invalpha3_c);

  sw2  = g(0)*g(0)*(3./5.)/(g(0)*g(0)*3./5. + g(1)*g(1));
  cw2  = 1. - sw2;


  ///////////////////////////////////////////////////////////////////////////
  // alpha_em(M_Z)_MSbar and Gmu
  ///////////////////////////////////////////////////////////////////////////
  double del_vb_sm, cw2pole;
  double del_vb_susy;
  //double a1, del_ve, del_vm, del_Zve, del_Zvm, del_Ze, del_Zm;
  double log_mu_mZ, log_md_mZ, log_me_mZ;
  MatrixXcd a_ncw(4,2), b_ncw(4,2);
  VectorXcd a_cev(2), b_cve(2), b_nvv(4), b_nee(4);

  // third family scalars contribution
  log_mu_mZ = 0.;
  log_md_mZ = 0.;
  log_me_mZ = 0.;
  /*
  for (int i=5; i<=6; i++) {
    log_mu_mZ = log_mu_mZ + log(mMu(i)/MEW);
    log_md_mZ = log_md_mZ + log(mMd(i)/MEW);
    log_me_mZ = log_me_mZ + log(mMe(i)/MEW);
  }
  */

  // alpha_3
  alpha_3 = ( -0.25 +
    (log_mu_mZ + log_md_mZ)/12. +
    //log(Mg/MEW) +
    log(mYu(2)/MEW)/3.)*g(2)*g(2)/(4.*pi*pi);
  alpha_3 = g(2)*g(2)/(4.*pi*( 1. - alpha_3));

  // alpha_em_DR
  alpha_em = sw2*g(1)*g(1)/(4.*pi);

  // Gmu : Nucl.Phys.B439,23
  Gmu = pi*alpha_em/(sqrt(2.)*(M_W2-PI_WW)*sw2);

  cw2pole = (M_W2 - PI_WW)/(M_Z2 - PI_ZZ);
  sw2pole = 1. - cw2pole;

  del_vb_sm = ( 6. + (7. - 5.*sw2pole + sw2*(3.*cw2pole/cw2 - 10.))
                     *log(cw2pole)/(2.*sw2pole) ) * alpha_em/(4.*pi*sw2);

  ///////////////////////////////////////////////////////////////////////////
  //delta_VB^SUSY temporarily removed for problems with C0 when masses equal!!!
  /*
  for (int i=1; i<=4; i++)
  {
    b_nvv(i) = - Nt(i,0)*g(0)*sqrt(3./10.) + Nt(i,1)*g(1)/sqrt(2.);
    b_nee(i) = - Nt(i,0)*g(0)*sqrt(3./10.) - Nt(i,1)*g(1)/sqrt(2.);
  }

  for (int i=1; i<=2; i++)
  {
    a_cev(i) = conj(V(i,0))*g(1);
    b_cve(i) = U(i,0)*g(1);
  }

  for (int i=1; i<=4; i++)
    for (int j=1; j<=2; j++)
    {
      a_ncw(i,j) = - g(1)*conj(Nt(i,1))*V(j,0)
                   + conj(Nt(i,3))*V(j,1)*g(1)/sqrt(2.);
      b_ncw(i,j) = - g(1)*Nt(i,1)*conj(U(j,0))
                   - Nt(i,2)*conj(U(j,1))*g(1)/sqrt(2.);
    }

  del_Zve = 0.;
  del_Zvm = 0.;
  del_Ze = 0.;
  del_Zm = 0.;
  del_ve = 0.;
  del_vm = 0.;
  a1 = 0.;
  p = 0.;

  //approximation!!!!!!!!!!!!!!!!!!!!!
  // selectron and sneutrino masses approximated by sqrt(ML2)
  // should be ok for first two generations, up to D-terms
  for (int i=1; i<=2; i++) {
    del_Zve = del_Zve + real(b_cve(i)*conj(b_cve(i)))
              *B1(p,mMch(i),sqrt(real(ML2(1,1))));
    del_Ze  = del_Ze  + real(a_cev(i)*conj(a_cev(i)))
              *B1(p,mMch(i),sqrt(real(ML2(1,1))));
    del_Zvm = del_Zvm + real(b_cve(i)*conj(b_cve(i)))
              *B1(p,mMch(i),sqrt(real(ML2(2,2))));
    del_Zm  = del_Zm  + real(a_cev(i)*conj(a_cev(i)))
              *B1(p,mMch(i),sqrt(real(ML2(2,2))));
  }

  for (int i=1; i<=4; i++) {
    del_Zve = del_Zve +
              real(b_nvv(i)*conj(b_nvv(i)))*B1(p,mMnt(i),sqrt(real(ML2(1,1))));
    del_Ze =  del_Ze  +
              real(b_nee(i)*conj(b_nee(i)))*B1(p,mMnt(i),sqrt(real(ML2(1,1))));
    del_Zvm = del_Zvm +
              real(b_nvv(i)*conj(b_nvv(i)))*B1(p,mMnt(i),sqrt(real(ML2(2,2))));
    del_Zm =  del_Zm  +
              real(b_nee(i)*conj(b_nee(i)))*B1(p,mMnt(i),sqrt(real(ML2(2,2))));
  }

  for (int i=1; i<=2; i++) {
    for (int j=1; j<=4; j++) {
      del_ve = del_ve + real( b_cve(i)*conj(b_nee(j)) * (
               - a_ncw(j,i) * mMch(i) * mMnt(j) *
               C0(sqrt(real(ML2(1,1))),mMch(i),mMnt(j))*sqrt(2.)/g(1) +
               b_ncw(j,i) * ( B0(p,mMch(i),mMnt(j)) + real(ML2(1,1)) *
               C0(sqrt(real(ML2(1,1))),mMch(i),mMnt(j))-0.5)/(sqrt(2.)*g(1)) ) )

               -  real( a_cev(i)*b_nvv(j) * (
               - b_ncw(j,i) * mMch(i) * mMnt(j) *
               C0(sqrt(real(ML2(1,1))),mMch(i),mMnt(j))*sqrt(2.)/g(1) +
               a_ncw(j,i) * ( B0(p,mMch(i),mMnt(j)) + real(ML2(1,1)) *
               C0(sqrt(real(ML2(1,1))),mMch(i),mMnt(j))-0.5)/(sqrt(2.)*g(1)) ) );

      del_vm = del_vm + real( b_cve(i)*conj(b_nee(j)) * (
               - a_ncw(j,i) * mMch(i) * mMnt(j) *
               C0(sqrt(real(ML2(2,2))),mMch(i),mMnt(j))*sqrt(2.)/g(1) +
               b_ncw(j,i) * ( B0(p,mMch(i),mMnt(j)) + real(ML2(2,2)) *
               C0(sqrt(real(ML2(2,2))),mMch(i),mMnt(j))-0.5)/(sqrt(2.)*g(1)) ) )

               -  real( a_cev(i)*b_nvv(j) * (
               - b_ncw(j,i) * mMch(i) * mMnt(j) *
               C0(sqrt(real(ML2(2,2))),mMch(i),mMnt(j))*sqrt(2.)/g(1) +
               a_ncw(j,i) * ( B0(p,mMch(i),mMnt(j)) + real(ML2(2,2)) *
               C0(sqrt(real(ML2(2,2))),mMch(i),mMnt(j))-0.5)/(sqrt(2.)*g(1)) ) );
    }
  }

  // neutralino contribution neglected at this point
  // because in this approximation C0 of two equal entry
  // causes problems
  for (int j=1; j<=4; j++) {
    del_ve = del_ve + 0.5* real( conj(b_nee(j))*b_nvv(j) ) * (B0(p,sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) + mMnt(j)*mMnt(j) *C0(mMnt(j),sqrt(real(ML2(1,1))),sqrt(real(ML2(1,1)))) + 0.5 );

    del_vm = del_vm + 0.5* real( conj(b_nee(j))*b_nvv(j) ) * (B0(p,sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) + mMnt(j)*mMnt(j) *C0(mMnt(j),sqrt(real(ML2(2,2))),sqrt(real(ML2(2,2)))) + 0.5 );
  }


  // only D0 terms nonzero in a1, D27 are zero for first two gens. degenerate!!!
  // second line of (C.20) is conjugate to the first line in this case !!!
  for (int i=1; i<=2; i++) {
    for (int j=1; j<=4; j++) {}
      a1 = a1 + real( a_cev(i) * conj(b_cve(i)) * b_nvv(j) * b_nee(j) ) * mMch(i)*mMnt(j) * D0(sqrt(real(ML2(1,1))),sqrt(real(ML2(2,2))),mMch(i),mMnt(j));
    }
  }

  del_vb_susy = ( - sw2*cw2*(M_Z2 - PI_ZZ)*a1/(2.*pi*alpha_em) + del_ve + del_vm + 0.5*( del_Zve + del_Zvm + del_Ze + del_Zm ) )/(16.*pi*pi);
  */

  del_vb_susy = 0.;

  double rho, del_rho, del_r;

  rho = cw2pole/cw2;
  del_rho = PI_ZZ/(rho*(M_Z2 - PI_ZZ)) - PI_WW/(M_W2 - PI_WW);

  del_r = PI_WW0/M_Wpole2 - PI_WW/M_Wpole2 + del_vb_sm + del_vb_susy;

  Gmu = Gmu/(1. - del_r);

  /////////////////////////////////////////////////////////////////////////////
  // alpha_em calculation : Pierce et. al. ////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  double factor;
  factor = (log_mu_mZ*4./9. + log_md_mZ/9. + log_me_mZ/3. +
    log(mMch(0)/MEW)*4./3. + log(mMch(1)/MEW)*4./3. +
    log(mHp/MEW)/3. +
    log(mYu(2)/MEW)*16./9. -
    log(M_Wpole2/(MEW*MEW))*7./2.) * alpha_em/(2.*pi);

  // alpha_em_MZ is corrected coupling in MS at MZ used for running of fermion masses below MZ!
  // different from Tomas! see sscor6.f - coupling
  // used for running nelow MZ is not corrected for top and W
  // and contains some constant term (from DR to MS???)!
  // Same as Pierce et. al.
  alpha_em_MZ = alpha_em /(1. - factor);

  // aplpa_em is the coupling at low energies should be compared with 1/137
  alpha_em = alpha_em /(1. - (- 0.0682 + factor));

  //cout<<"alpha_em with 0.0689 !!!"<<endl;
  //cout<<"mYu(2) = "<<mYu(2)<<endl;
  //cout<<"alpha_em with mYu(2)!!!"<<endl;

  return Pe;
}
