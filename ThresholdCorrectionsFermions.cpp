#include "ThresholdCorrectionsFermions.hpp"
#include "PassarinoVeltmanFunctions.h"
//-----------------------------------------------------------------------
// Calculating threshold corrections to quarks and leptons at Z scale
// output is written in Yukava matrices, so the Yu contains original mass
// eigenstates + corrections
//-----------------------------------------------------------------------

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
        double M_Z, double M_W, double sw2, double mA, double mHp, double mH,
        double mh, double c2a, double s2a)
{
    
  const double pi = 3.141592653589793;
  const double sixteenpi2 = 16.* pi * pi;

  double pt = 174.3;
  double p0 = 0.;
 
// gluon correction to top

  double dm_top_g, vu, vd, Q;

//*****************************************************************************8
// variables for gluino corrections
  double mg, fac_gino;

  MatrixXcd
                dL_d_gino(3,3), dR_d_gino(3,3), dm_d_gino(3,3), dm_d_ginot(3,3),
                dL_u_gino(3,3), dR_u_gino(3,3), dm_u_gino(3,3), dm_u_ginot(3,3);

//*****************************************************************************8
// variables for chargino corrections

  std::complex<double> ChdR[2][6][3], ChdL[2][6][3], ChuR[2][6][3], ChuL[2][6][3],
                  CheR[2][3][3], CheL[2][3][3];

  MatrixXcd
                      dL_d_ch(3,3), dR_d_ch(3,3), dm_d_ch(3,3), dm_d_cht(3,3),
                      dL_u_ch(3,3), dR_u_ch(3,3), dm_u_ch(3,3), dm_u_cht(3,3),
                      dL_e_ch(3,3), dR_e_ch(3,3), dm_e_ch(3,3), dm_e_cht(3,3);

//*****************************************************************************8
// variables for neutralino corrections
        
  std::complex<double> NtdR[4][6][3], NtdL[4][6][3], NtuR[4][6][3], NtuL[4][6][3],
                  NteR[4][6][3], NteL[4][6][3];

  MatrixXcd
                     dL_d_nt(3,3), dR_d_nt(3,3), dm_d_nt(3,3), dm_d_ntt(3,3),
                     dL_u_nt(3,3), dR_u_nt(3,3), dm_u_nt(3,3), dm_u_ntt(3,3),
                     dL_e_nt(3,3), dR_e_nt(3,3), dm_e_nt(3,3), dm_e_ntt(3,3);

//*****************************************************************************8
// variables for Charged Higgs corrections

  MatrixXcd
                      dL_d_Hp(3,3), dR_d_Hp(3,3), dm_d_Hp(3,3), dm_d_Hpt(3,3),
                      dL_u_Hp(3,3), dR_u_Hp(3,3), dm_u_Hp(3,3), dm_u_Hpt(3,3),
                      dL_e_Hp(3,3), dR_e_Hp(3,3), dm_e_Hp(3,3), dm_e_Hpt(3,3);

//*****************************************************************************8
// variables for Neutral Higgs corrections

  MatrixXcd
                      dL_d_H0(3,3), dR_d_H0(3,3), dm_d_H0(3,3), dm_d_H0t(3,3),
                      dL_u_H0(3,3), dR_u_H0(3,3), dm_u_H0(3,3), dm_u_H0t(3,3),
                      dL_e_H0(3,3), dR_e_H0(3,3), dm_e_H0(3,3), dm_e_H0t(3,3);

//*****************************************************************************8
// variables for W corrections

  MatrixXcd
                      dL_d_W(3,3), dR_d_W(3,3), dm_d_W(3,3), dm_d_Wt(3,3),
                      dL_u_W(3,3), dR_u_W(3,3), dm_u_W(3,3), dm_u_Wt(3,3),
                      dL_e_W(3,3), dR_e_W(3,3), dm_e_W(3,3), dm_e_Wt(3,3);

//*****************************************************************************8
// variables for Z corrections

  MatrixXcd
                      dL_d_Z(3,3), dR_d_Z(3,3), dm_d_Z(3,3), dm_d_Zt(3,3),
                      dL_u_Z(3,3), dR_u_Z(3,3), dm_u_Z(3,3), dm_u_Zt(3,3),
                      dL_e_Z(3,3), dR_e_Z(3,3), dm_e_Z(3,3), dm_e_Zt(3,3);

//*****************************************************************************8
// Higgs, Z, W masses and mixing angles, running DRbar at M_Z

  double cb, sb, cb2, sb2, sa2, ca2;

  cb  = 1./sqrt(1.+tanb*tanb);
  sb  = tanb*cb;
  cb2 = cb*cb;
  sb2 = 1. - cb2;

  sa2 = (1. - c2a)/2.;
  ca2 = 1. - sa2;

  vu = v / sqrt(2. + 2. /(tanb*tanb));
  vd = v / sqrt(2. + 2.* tanb * tanb);
  Q = 91.188;

for (int i=0; i<3; i++)
{
  mYu(i) = mYu(i)*(1./vu);
  mYd(i) = mYd(i)*(1./vd);
  mYe(i) = mYe(i)*(1./vd);
}

//*****************************************************************************8
// gluino corrections

  mg = std::real(M(2));

//  mg = sqrt(real(M[2])*real(M[2]) + imag(M[2])*imag(M[2]));
//  mg = 789.0928;

  fac_gino = g[2] * g[2] * 8. /(3. * sixteenpi2);

  for (int i = 0; i <= 2; i++)
  {
    for (int j = 0; j <= 2; j++)
    {
      dL_d_gino(i, j) = dcomplex(0., 0.);
      dR_d_gino(i, j) = dcomplex(0., 0.);
      dm_d_gino(i, j) = dcomplex(0., 0.);

      dL_u_gino(i, j) = dcomplex(0., 0.);
      dR_u_gino(i, j) = dcomplex(0., 0.);
      dm_u_gino(i, j) = dcomplex(0., 0.);

      for (int l = 0; l <= 5; l++)
      {
        dL_d_gino(i, j) = dL_d_gino(i, j) +
                          fac_gino*conj(GdL(l, i))*GdL(l, j)*B1(p0,mg,mMd[l]);
        dR_d_gino(i, j) = dR_d_gino(i, j) +
                          fac_gino*conj(GdR(l, i))*GdR(l, j)*B1(p0,mg,mMd[l]);
        dm_d_gino(i, j) = dm_d_gino(i, j) - mg *
                          fac_gino*conj(GdL(l, i))*GdR(l, j)*B0(p0,mg,mMd[l]);

        dL_u_gino(i, j) = dL_u_gino(i, j) +
                          fac_gino*conj(GuL(l, i))*GuL(l, j)*B1(pt,mg,mMu[l]);
        dR_u_gino(i, j) = dR_u_gino(i, j) +
                          fac_gino*conj(GuR(l, i))*GuR(l, j)*B1(pt,mg,mMu[l]);
        dm_u_gino(i, j) = dm_u_gino(i, j) - mg*
                          fac_gino*conj(GuL(l, i))*GuR(l, j)*B0(pt,mg,mMu[l]);
      }

      dm_d_ginot(i, j) = .5*vd*(dL_d_gino(i, j)*mYd(j) + mYd[i]*dR_d_gino(i, j))
                         - dm_d_gino(i, j);

      dm_u_ginot(i, j) = .5*vu*(dL_u_gino(i, j)*mYu(j) + mYu[i]*dR_u_gino(i, j))
                         - dm_u_gino(i, j);
    }
  }
/*
cout<<"dL_u_gino: "<<dL_u_gino<<endl;
cout<<"dR_u_gino: "<<dR_u_gino<<endl;
cout<<"dm_u_gino: "<<dm_u_gino<<endl;


cout<<"dL_d_gino: "<<dL_d_gino<<endl;
cout<<"dR_d_gino: "<<dR_d_gino<<endl;
cout<<"dm_d_gino: "<<dm_d_gino<<endl;
*/

/*
// gluon correction to top quark mass

  dm_top_g = vu*mYu[2] * .5*fac_gino * (6.*log(Q/(mYu[2]*vu)) + 5.);

// cout<<"EtaTop: "<< (mYu(3) + dm_top_g)/mYu(3)<<endl;
*/

//*****************************************************************************8
// chargino corrections

  for (int i = 0; i <= 2; i++)
  {
    for (int l = 0; l <= 5; l++)
    {
      for (int A = 0; A <= 1; A++)
      {
        ChdL[A][l][i] = dcomplex(0., 0.);
        ChdR[A][l][i] = dcomplex(0., 0.);
        ChuL[A][l][i] = dcomplex(0., 0.);
        ChuR[A][l][i] = dcomplex(0., 0.);

        for (int j = 0; j <= 2; j++)
        {
          ChdL[A][l][i] = ChdL[A][l][i] + ( -g[1]*conj(V(A, 0))*GuL(l, j) +
                          mYu[j]*conj(V(A, 1))*GuR(l, j) ) * VCKM(j, i);

          ChdR[A][l][i] = ChdR[A][l][i] +
                          mYd[i]*U(A, 1)*GuL(l, j)*VCKM(j, i);

          ChuL[A][l][i] = ChuL[A][l][i] + (-g[1]*conj(U(A, 0))*GdL(l, j) +
                          mYd[j]*conj(U(A, 1))*GdR(l, j))*conj(VCKM(i, j));

          ChuR[A][l][i] = ChuR[A][l][i] +
                          mYu[i]*V(A, 1)*GdL(l, j)*conj(VCKM(i, j));

        }
      }
    }
  }

  for (int i = 0; i <= 2; i++)
  {
    for (int l = 0; l <= 2; l++)
    {
      for (int A = 0; A <= 1; A++)
      {
        CheL[A][l][i] = - g[1]*conj(V(A, 0))*GvL(l, i);
        CheR[A][l][i] = mYe[i]*U(A, 1)*GvL(l, i);

      }
    }
  }

  for (int i = 0; i <= 2; i++)
  {
    for (int j = 0; j <= 2; j++)
    {
      dL_d_ch(i, j) = dcomplex(0., 0.);
      dR_d_ch(i, j) = dcomplex(0., 0.);
      dm_d_ch(i, j) = dcomplex(0., 0.);

      dL_u_ch(i, j) = dcomplex(0., 0.);
      dR_u_ch(i, j) = dcomplex(0., 0.);
      dm_u_ch(i, j) = dcomplex(0., 0.);
             
      dL_e_ch(i, j) = dcomplex(0., 0.);
      dR_e_ch(i, j) = dcomplex(0., 0.);
      dm_e_ch(i, j) = dcomplex(0., 0.);
     
      for (int l = 0; l <= 5; l++)
      {
        for (int A = 0; A <= 1; A++)
        {
          dL_d_ch(i, j) = dL_d_ch(i, j) + conj(ChdL[A][l][i])*ChdL[A][l][j]
                          * B1(p0,mMch[A],mMu[l])/sixteenpi2;
          dR_d_ch(i, j) = dR_d_ch(i, j) + conj(ChdR[A][l][i])*ChdR[A][l][j]
                          * B1(p0,mMch[A],mMu[l])/sixteenpi2;
          dm_d_ch(i, j) = dm_d_ch(i, j) + conj(ChdL[A][l][i])*ChdR[A][l][j]
                          * B0(p0,mMch[A],mMu[l])*mMch[A]/sixteenpi2;

          dL_u_ch(i, j) = dL_u_ch(i, j) + conj(ChuL[A][l][i])*ChuL[A][l][j]
                          * B1(pt,mMch[A],mMd[l])/sixteenpi2;
          dR_u_ch(i, j) = dR_u_ch(i, j) + conj(ChuR[A][l][i])*ChuR[A][l][j]
                          * B1(pt,mMch[A],mMd[l])/sixteenpi2;
          dm_u_ch(i, j) = dm_u_ch(i, j) + conj(ChuL[A][l][i])*ChuR[A][l][j]
                          * B0(pt,mMch[A],mMd[l])*mMch[A]/sixteenpi2;
        }
      }

      for (int l = 0; l <= 2; l++)
      {
        for (int A = 0; A <= 1; A++)
        {
          dL_e_ch(i, j) = dL_e_ch(i, j) + conj(CheL[A][l][i])*CheL[A][l][j]
                          * B1(p0,mMch[A],0)/sixteenpi2;
          dR_e_ch(i, j) = dR_e_ch(i, j) + conj(CheR[A][l][i])*CheR[A][l][j]
                          * B1(p0,mMch[A],0)/sixteenpi2;
          dm_e_ch(i, j) = dm_e_ch(i, j) + conj(CheL[A][l][i])*CheR[A][l][j]
                          * B0(p0,mMch[A],0)*mMch[A]/sixteenpi2;
        }
      }

      dm_d_cht(i, j) = .5*vd*(dL_d_ch(i, j)*mYd[j] + mYd[i]*dR_d_ch(i, j))
                       - dm_d_ch(i, j);
      dm_u_cht(i, j) = .5*vu*(dL_u_ch(i, j)*mYu[j] + mYu[i]*dR_u_ch(i, j))
                       - dm_u_ch(i, j);
      dm_e_cht(i, j) = .5*vd*(dL_e_ch(i, j)*mYe[j] + mYe[i]*dR_e_ch(i, j))
                       - dm_e_ch(i, j);
    }
  }
/*
cout<<"dL_u_ch: "<<dL_u_ch<<endl;
cout<<"dR_u_ch: "<<dR_u_ch<<endl;
cout<<"dm_u_ch: "<<dm_u_ch<<endl;

cout<<"dL_d_ch: "<<dL_d_ch<<endl;
cout<<"dR_d_ch: "<<dR_d_ch<<endl;
cout<<"dm_d_ch: "<<dm_d_ch<<endl;

cout<<"dL_e_ch: "<<dL_e_ch<<endl;
cout<<"dR_e_ch: "<<dR_e_ch<<endl;
cout<<"dm_e_ch: "<<dm_e_ch<<endl;
*/

//*****************************************************************************8
// neutralino corrections

  double gp = g[0]*sqrt(3./5.);

  for (int i = 0; i <= 2; i++)
  {
    for (int l = 0; l <= 5; l++)
    {
      for (int A = 0; A <= 3; A++)
      {
        NtdL[A][l][i] = ( g[1]*conj(Nt(A, 1)) - gp*conj(Nt(A, 0))/3. )
                        * GdL(l, i)/sqrt(2.) - mYd[i]*conj(Nt(A, 2))*GdR(l, i);
        NtdR[A][l][i] = (-gp*Nt(A, 0)*2./3.)*GdR(l, i)/sqrt(2.)
                        - mYd[i]*Nt(A, 2)*GdL(l, i);

        NtuL[A][l][i] = (-g[1]*conj(Nt(A, 1)) - gp*conj(Nt(A, 0))/3. )
                        * GuL(l, i)/sqrt(2.) - mYu[i]*conj(Nt(A, 3))*GuR(l, i);
        NtuR[A][l][i] = (gp*Nt(A, 0)*4./3.)*GuR(l, i)/sqrt(2.)
                        - mYu[i]*Nt(A, 3)*GuL(l, i);

        NteL[A][l][i] = ( g[1]*conj(Nt(A, 1)) + gp*conj(Nt(A, 0)) )
                        * GeL(l, i)/sqrt(2.) - mYe[i]*conj(Nt(A, 2))*GeR(l, i);
        NteR[A][l][i] = (-gp*Nt(A, 0)*2.)*GeR(l, i)/sqrt(2.)
                        - mYe[i]*Nt(A, 2)*GeL(l, i);
      }
    }
  }

  for (int i = 0; i <= 2; i++)
  {
    for (int j = 0; j <= 2; j++)
    {
      dL_d_nt(i, j) = dcomplex(0., 0.);
      dR_d_nt(i, j) = dcomplex(0., 0.);
      dm_d_nt(i, j) = dcomplex(0., 0.);

      dL_u_nt(i, j) = dcomplex(0., 0.);
      dR_u_nt(i, j) = dcomplex(0., 0.);
      dm_u_nt(i, j) = dcomplex(0., 0.);

      dL_e_nt(i, j) = dcomplex(0., 0.);
      dR_e_nt(i, j) = dcomplex(0., 0.);
      dm_e_nt(i, j) = dcomplex(0., 0.);
 
      for (int l = 0; l <= 5; l++)
      {
        for (int A = 0; A <= 3; A++)
        {
          dL_d_nt(i, j) = dL_d_nt(i, j) + conj(NtdL[A][l][i])*NtdL[A][l][j]
                          * B1(p0,mMnt[A],mMd[l])/sixteenpi2;
          dR_d_nt(i, j) = dR_d_nt(i, j) + conj(NtdR[A][l][i])*NtdR[A][l][j]
                          * B1(p0,mMnt[A],mMd[l])/sixteenpi2;
          dm_d_nt(i, j) = dm_d_nt(i, j) + conj(NtdL[A][l][i])*NtdR[A][l][j]
                          * B0(p0,mMnt[A],mMd[l])*mMnt[A]/sixteenpi2;

          dL_u_nt(i, j) = dL_u_nt(i, j) + conj(NtuL[A][l][i])*NtuL[A][l][j]
                          * B1(pt,mMnt[A],mMu[l])/sixteenpi2;
          dR_u_nt(i, j) = dR_u_nt(i, j) + conj(NtuR[A][l][i])*NtuR[A][l][j]
                          * B1(pt,mMnt[A],mMu[l])/sixteenpi2;
          dm_u_nt(i, j) = dm_u_nt(i, j) + conj(NtuL[A][l][i])*NtuR[A][l][j]
                          * B0(pt,mMnt[A],mMu[l])*mMnt[A]/sixteenpi2;
  
          dL_e_nt(i, j) = dL_e_nt(i, j) + conj(NteL[A][l][i])*NteL[A][l][j]
                          * B1(p0,mMnt[A],mMe[l])/sixteenpi2;
          dR_e_nt(i, j) = dR_e_nt(i, j) + conj(NteR[A][l][i])*NteR[A][l][j]
                          * B1(p0,mMnt[A],mMe[l])/sixteenpi2;
          dm_e_nt(i, j) = dm_e_nt(i, j) + conj(NteL[A][l][i])*NteR[A][l][j]
                          * B0(p0,mMnt[A],mMe[l])*mMnt[A]/sixteenpi2;
        }
      }

      dm_d_ntt(i, j) = .5*vd*(dL_d_nt(i, j)*mYd[j] + mYd[i]*dR_d_nt(i, j))
                       - dm_d_nt(i, j);
      dm_u_ntt(i, j) = .5*vu*(dL_u_nt(i, j)*mYu[j] + mYu[i]*dR_u_nt(i, j))
                       - dm_u_nt(i, j);
      dm_e_ntt(i, j) = .5*vd*(dL_e_nt(i, j)*mYe[j] + mYe[i]*dR_e_nt(i, j))
                       - dm_e_nt(i, j);

    }
  }


//*****************************************************************************8
// Charged Higgs corrections

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      dL_d_Hp(i,j) = dcomplex(0., 0.);
      dR_d_Hp(i,j) = dcomplex(0., 0.);
      dm_d_Hp(i,j) = dcomplex(0., 0.);

      dL_u_Hp(i,j) = dcomplex(0., 0.);
      dR_u_Hp(i,j) = dcomplex(0., 0.);
      dm_u_Hp(i,j) = dcomplex(0., 0.);

      dL_e_Hp(i,j) = dcomplex(0., 0.);
      dR_e_Hp(i,j) = dcomplex(0., 0.);
      dm_e_Hp(i,j) = dcomplex(0., 0.);

      for (int k = 0; k < 3; k++)
      {
        dL_d_Hp(i,j) = dL_d_Hp(i,j) + conj(VCKM(k,i))*mYu(k)*mYu(k) *
                        ( B1(p0,mYu(k)*vu,mHp)*cb2 + B1(p0,mYu(k)*vu,M_W)*sb2 )
                        * VCKM(k,j)/sixteenpi2;

        dR_d_Hp(i,j) = dR_d_Hp(i,j) + mYd(i)*conj(VCKM(k,i))*
                        ( B1(p0,mYu(k)*vu,mHp)*sb2 + B1(p0,mYu(k)*vu,M_W)*cb2 )
                        * VCKM(k,j)*mYd(j)/sixteenpi2;

        dm_d_Hp(i,j) = dm_d_Hp(i,j) + conj(VCKM(k,i))*mYu(k)*mYu(k)*vu*
                       ( B0(p0,mYu(k)*vu,mHp) - B0(p0,mYu(k)*vu,M_W) )
                       *cb*sb* VCKM(k,j)*mYd(j)/sixteenpi2;


        dL_u_Hp(i,j) = dL_u_Hp(i,j) + VCKM(i,k)*mYd(k)*mYd(k) *
                        ( B1(pt,mYd(k)*vd,mHp)*sb2 + B1(pt,mYd(k)*vd,M_W)*cb2 )
                        * conj(VCKM(j,k))/sixteenpi2;

        dR_u_Hp(i,j) = dR_u_Hp(i,j) + mYu(i)*VCKM(i,k)*
                        ( B1(pt,mYd(k)*vd,mHp)*cb2 + B1(pt,mYd(k)*vd,M_W)*sb2 )
                        * conj(VCKM(j,k))*mYu(j)/sixteenpi2;

        dm_u_Hp(i,j) = dm_u_Hp(i,j) + VCKM(i,k)*mYd(k)*mYd(k)*vd*
                       ( B0(pt,mYd(k)*vd,mHp) - B0(pt,mYd(k)*vd,M_W) )
                       *cb*sb* conj(VCKM(j,k))*mYu(j)/sixteenpi2;

      }

      dm_d_Hpt(i,j) = .5*vd*(dL_d_Hp(i,j)*mYd(j) + mYd(i)*dR_d_Hp(i,j))
                       - dm_d_Hp(i,j);
      dm_u_Hpt(i,j) = .5*vu*(dL_u_Hp(i,j)*mYu(j) + mYu(i)*dR_u_Hp(i,j))
                       - dm_u_Hp(i,j);

      dm_e_Hpt(i,j) = dcomplex(0.,0.);
    }

    dR_e_Hp(i,i) = mYe(i)*mYe(i)*
                   ( B1(p0,p0,mHp)*sb2 + B1(p0,p0,M_W)*cb2 ) /sixteenpi2;

    dm_e_Hpt(i,i) = .5*vd * mYe(i)*dR_e_Hp(i,i);
  }

//*****************************************************************************8
// Neutral Higgs corrections

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      dL_d_H0(i,j) = std::complex<double> (0., 0.);
      dR_d_H0(i,j) = std::complex<double> (0., 0.);
      dm_d_H0(i,j) = std::complex<double> (0., 0.);

      dL_u_H0(i,j) = std::complex<double> (0., 0.);
      dR_u_H0(i,j) = std::complex<double> (0., 0.);
      dm_u_H0(i,j) = std::complex<double> (0., 0.);

      dL_e_H0(i,j) = std::complex<double> (0., 0.);
      dR_e_H0(i,j) = std::complex<double> (0., 0.);
      dm_e_H0(i,j) = std::complex<double> (0., 0.);

    }

    dL_d_H0(i,i) = 0.5*mYd(i)*mYd(i) *
                   ( B1(p0,mYd(i)*vd,mH)*ca2 + B1(p0,mYd(i)*vd,mh)*sa2 +
                     B1(p0,mYd(i)*vd,mA)*sb2 + B1(p0,mYd(i)*vd,M_Z)*cb2 )
                   /sixteenpi2;

    dR_d_H0(i,i) = dL_d_H0(i,i);

    dm_d_H0(i,i) = 0.5*mYd(i)*mYd(i)*mYd(i)*vd *
                   ( B0(p0,mYd(i)*vd,mH)*ca2 + B0(p0,mYd(i)*vd,mh)*sa2 -
                     B0(p0,mYd(i)*vd,mA)*sb2 - B0(p0,mYd(i)*vd,M_Z)*cb2 )
                   /sixteenpi2;


    dL_u_H0(i,i) = 0.5*mYu(i)*mYu(i) *
                   ( B1(pt,mYu(i)*vu,mH)*sa2 + B1(pt,mYu(i)*vu,mh)*ca2 +
                     B1(pt,mYu(i)*vu,mA)*cb2 + B1(pt,mYu(i)*vu,M_Z)*sb2 )
                   /sixteenpi2;

    dR_u_H0(i,i) = dL_u_H0(i,i);

    dm_u_H0(i,i) = 0.5*mYu(i)*mYu(i)*mYu(i)*vu *
                   ( B0(pt,mYu(i)*vu,mH)*sa2 + B0(pt,mYu(i)*vu,mh)*ca2 -
                     B0(pt,mYu(i)*vu,mA)*cb2 - B0(pt,mYu(i)*vu,M_Z)*sb2 )
                   /sixteenpi2;


    dL_e_H0(i,i) = 0.5*mYe(i)*mYe(i) *
                   ( B1(p0,mYe(i)*vd,mH)*ca2 + B1(p0,mYe(i)*vd,mh)*sa2 +
                     B1(p0,mYe(i)*vd,mA)*sb2 + B1(p0,mYe(i)*vd,M_Z)*cb2 )
                   /sixteenpi2;

    dR_e_H0(i,i) = dL_e_H0(i,i);

    dm_e_H0(i,i) = 0.5*mYe(i)*mYe(i)*mYe(i)*vd *
                   ( B0(p0,mYe(i)*vd,mH)*ca2 + B0(p0,mYe(i)*vd,mh)*sa2 -
                     B0(p0,mYe(i)*vd,mA)*sb2 - B0(p0,mYe(i)*vd,M_Z)*cb2 )
                   /sixteenpi2;

  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      
      dm_d_H0t(i,j) = .5*vd*(dL_d_H0(i,j)*mYd(j) + mYd(i)*dR_d_H0(i,j))
                       - dm_d_H0(i,j);

      dm_u_H0t(i,j) = .5*vu*(dL_u_H0(i,j)*mYu(j) + mYu(i)*dR_u_H0(i,j))
                       - dm_u_H0(i,j);

      dm_e_H0t(i,j) = .5*vd*(dL_e_H0(i,j)*mYe(j) + mYe(i)*dR_e_H0(i,j))
                       - dm_e_H0(i,j);
    }
  }


//*****************************************************************************8
// W corrections

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      dL_d_W(i,j) = std::complex<double> (0., 0.);
      dR_d_W(i,j) = std::complex<double> (0., 0.);
      dm_d_W(i,j) = std::complex<double> (0., 0.);

      dL_u_W(i,j) = std::complex<double> (0., 0.);
      dR_u_W(i,j) = std::complex<double> (0., 0.);
      dm_u_W(i,j) = std::complex<double> (0., 0.);

      dL_e_W(i,j) = std::complex<double> (0., 0.);
      dR_e_W(i,j) = std::complex<double> (0., 0.);
      dm_e_W(i,j) = std::complex<double> (0., 0.);

      for (int k = 0; k < 3; k++)
      {
        dL_d_W(i,j) = dL_d_W(i,j) + g(1)*g(1) * conj(VCKM(k,i)) *
                      B1(p0,mYu(k)*vu,M_W) * VCKM(k,j) /sixteenpi2;

        dL_u_W(i,j) = dL_u_W(i,j) + g(1)*g(1) * VCKM(i,k) *
                      B1(pt,mYd(k)*vd,M_W) * conj(VCKM(j,k)) /sixteenpi2;

      }
    }

    dL_e_W(i,i) = dL_e_W(i,i) + g(1)*g(1) * B1(p0,p0,M_W) /sixteenpi2;

  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      dm_d_Wt(i,j) = .5*vd*(dL_d_W(i,j)*mYd(j) + mYd(i)*dR_d_W(i,j))
                       - dm_d_W(i,j);

      dm_u_Wt(i,j) = .5*vu*(dL_u_W(i,j)*mYu(j) + mYu(i)*dR_u_W(i,j))
                       - dm_u_W(i,j);

      dm_e_Wt(i,j) = .5*vd*(dL_e_W(i,j)*mYe(j) + mYe(i)*dR_e_W(i,j))
                       - dm_e_W(i,j);
    }
  }

//*****************************************************************************8
// Z corrections

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      dL_d_Z(i,j) = std::complex<double> (0., 0.);
      dR_d_Z(i,j) = std::complex<double> (0., 0.);
      dm_d_Z(i,j) = std::complex<double> (0., 0.);

      dL_u_Z(i,j) = std::complex<double> (0., 0.);
      dR_u_Z(i,j) = std::complex<double> (0., 0.);
      dm_u_Z(i,j) = std::complex<double> (0., 0.);

      dL_e_Z(i,j) = std::complex<double> (0., 0.);
      dR_e_Z(i,j) = std::complex<double> (0., 0.);
      dm_e_Z(i,j) = std::complex<double> (0., 0.);
    }

    dL_d_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*(-0.5 + sw2/3.)*(-0.5 + sw2/3.)*B1(p0,mYd(i)*vd,M_Z)
                   /sixteenpi2;

    dR_d_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*(sw2/3.)*(sw2/3.)*B1(p0,mYd(i)*vd,M_Z)
                   /sixteenpi2;

    dm_d_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) * mYd(i)*vd *
                   4.*(-0.5 + sw2/3.)*(-sw2/3.)*B0(p0,mYd(i)*vd,M_Z)
                   /sixteenpi2;


    dL_u_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*(0.5 - sw2*2./3.)*(0.5 - sw2*2./3.)*B1(pt,mYu(i)*vu,M_Z)
                   /sixteenpi2;

    dR_u_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*(sw2*2./3.)*(sw2*2./3.)*B1(pt,mYu(i)*vu,M_Z)
                   /sixteenpi2;

    dm_u_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) * mYu(i)*vu *
                   4.*(0.5 - sw2*2./3.)*(sw2*2./3.)*B0(pt,mYu(i)*vu,M_Z)
                   /sixteenpi2;


    dL_e_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*(-0.5 + sw2)*(-0.5 + sw2)*B1(p0,mYe(i)*vd,M_Z)
                   /sixteenpi2;

    dR_e_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) *
                   2.*sw2*sw2*B1(p0,mYe(i)*vd,M_Z)
                   /sixteenpi2;

    dm_e_Z(i,i) =  (g(0)*g(0)*3./5 + g(1)*g(1)) * mYe(i)*vd *
                   4.*(-0.5 + sw2)*(-sw2)*B0(p0,mYe(i)*vd,M_Z)
                   /sixteenpi2;
  }

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {

      dm_d_Zt(i,j) = .5*vd*(dL_d_Z(i,j)*mYd(j) + mYd(i)*dR_d_Z(i,j))
                       - dm_d_Z(i,j);

      dm_u_Zt(i,j) = .5*vu*(dL_u_Z(i,j)*mYu(j) + mYu(i)*dR_u_Z(i,j))
                       - dm_u_Z(i,j);

      dm_e_Zt(i,j) = .5*vd*(dL_e_Z(i,j)*mYe(j) + mYe(i)*dR_e_Z(i,j))
                       - dm_e_Z(i,j);
    }
  }

/*
cout<<"vu = "<<vu<<endl;
cout<<"vd = "<<vd<<endl;
cout<<"mYu = "<<mYu<<endl;
cout<<"Yu = "<<Yu<<endl;
*/
//*****************************************************************************8
// after the last correction is calculated all corrections are summed up

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      Yd(i,j) =  dm_d_ginot(i,j) + dm_d_cht(i,j) + dm_d_ntt(i,j) + dm_d_Hpt(i,j)
               + dm_d_H0t(i,j)   + dm_d_Wt(i,j)  + dm_d_Zt(i,j);

      Yu(i,j) =  dm_u_ginot(i,j) + dm_u_cht(i,j) + dm_u_ntt(i,j) + dm_u_Hpt(i,j)
               + dm_u_H0t(i,j)   + dm_u_Wt(i,j)  + dm_u_Zt(i,j);

      Ye(i,j) =                    dm_e_cht(i,j) + dm_e_ntt(i,j) + dm_e_Hpt(i,j)
               + dm_e_H0t(i,j)   + dm_e_Wt(i,j)  + dm_e_Zt(i,j);
    }

// and finally new quark and lepton matrices are formed

    Yd(i,i) = Yd(i,i) + vd*mYd(i);
    Yu(i,i) = Yu(i,i) + vu*mYu(i);
    Ye(i,i) = Ye(i,i) + vd*mYe(i);
  }

//*****************************************************************************8
// in my notation W= Q lam u, or equivalently L = - barQ lam^* u
// corrections calculated here are actually corrections to lam^*
// so we conjugate matrices here

  Yu = Yu.conjugate().eval();
  Yd = Yd.conjugate().eval();
  Ye = Ye.conjugate().eval();

}
