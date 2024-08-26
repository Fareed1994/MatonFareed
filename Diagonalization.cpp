#include "Diagonalization.hpp"
#include "PassarinoVeltmanFunctions.h"
#include "ThresholdCorrectionsFermions.hpp"
//---------------------------------------------------------------------
// diagonalization of (s)quarks, (s)leptons, charginos and neutralinos
//---------------------------------------------------------------------

void order(int n, VectorXd& m, MatrixXcd& A)
{
    // for the time being it works for n=3

    VectorXd mzal(3);
    MatrixXcd Azal(3, 3);

    mzal = m;
    Azal = A;

    if (mzal(0) >= mzal(1) && mzal(0) >= mzal(2))
    {
        if (mzal(1) >= mzal(2))
        {
            m(2) = mzal(0);
            m(0) = mzal(2);
            for (int i = 0; i < n; i++)
            {
                A(i, 2) = Azal(i, 0);
                A(i, 0) = Azal(i, 2);
            }
        }
        else
        {
            m(2) = mzal(0);
            m(1) = mzal(2);
            m(0) = mzal(1);
            for (int i = 0; i < n; i++)
            {
                A(i, 2) = Azal(i, 0);
                A(i, 1) = Azal(i, 2);
                A(i, 0) = Azal(i, 1);
            }
        }
    }
    else if (mzal(1) >= mzal(0) && mzal(1) >= mzal(2))
    {
        if (mzal(2) >= mzal(0))
        {
            m(2) = mzal(1);
            m(1) = mzal(2);
            for (int i = 0; i < n; i++)
            {
                A(i, 2) = Azal(i, 1);
                A(i, 1) = Azal(i, 2);
            }
        }
        else
        {
            m(2) = mzal(1);
            m(1) = mzal(0);
            m(0) = mzal(2);
            for (int i = 0; i < n; i++)
            {
                A(i, 2) = Azal(i, 1);
                A(i, 1) = Azal(i, 0);
                A(i, 0) = Azal(i, 2);
            }
        }
    }
    else
    {
        m(0) = mzal(1);
        m(1) = mzal(0);
        for (int i = 0; i < n; i++)
        {
            A(i, 0) = Azal(i, 1);
            A(i, 1) = Azal(i, 0);
        }
    }

    return;
}


// ###########################################################################
// reduction of a std::real symmetric matrix to tridiagonal form
// Input:  a[n][n] is a std::real symmetric matrix
// Output: a[n][n] is orthogonal matrix Q effecting the transformation, it is
//                 an input for the subroutine tqli which finish diag.
//         d[n] contains diagonal elements of the tridiagonal matrix,
//         e[n] contains subdiagonal elements, with e(1) = 0.
// ###########################################################################
void tred2(MatrixXd& a, int n, VectorXd& d, VectorXd& e)
{
  int l, k, j, i;
  double scale, hh, h, g, f;

  for (i = n - 1; i >= 1; i--)
  {
    l = i - 1;
    h = scale = 0.0;
    if (l > 0)
    {
      for (k = 0; k <= l; k++)
        scale += fabs(a(i, k));
      if (scale == 0.0)
        e(i) = a(i, l);
      else
      {
        for (k = 0; k <= l; k++)
        {
          a(i, k) /= scale;
          h += a(i, k) * a(i, k);
        }
        f = a(i, l);
        g = (f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
        e(i) = scale * g;
        h -= f * g;
        a(i, l) = f - g;
        f = 0.0;
        for (j = 0; j <= l; j++)
        {
          a(j, i) = a(i, j) / h;
          g = 0.0;
          for (k = 0; k <= j; k++)
            g += a(j, k) * a(i, k);
          for (k = j + 1; k <= l; k++)
            g += a(k, j) * a(i, k);
          e(j) = g / h;
          f += e(j) * a(i, j);
        }
        hh = f / (h + h);
        for (j = 0; j <= l; j++)
        {
          f = a(i, j);
          e(j) = g = e(j) - hh * f;
          for (k = 0; k <= j; k++)
            a(j, k) -= (f * e(k) + g * a(i, k));
        }
      }
    }
    else
      e(i) = a(i, l);
    d(i) = h;
  }
  d(0) = 0.0;
  e(0) = 0.0;

  for (i = 0; i < n; i++)
  {
    l = i - 1;
    if (d(i))
    {
      for (j = 0; j <= l; j++)
      {
        g = 0.0;
        for (k = 0; k <= l; k++)
          g += a(i, k) * a(k, j);
        for (k = 0; k <= l; k++)
          a(k, j) -= g * a(k, i);
      }
    }
    d(i) = a(i, i);
    a(i, i) = 1.0;
    for (j = 0; j <= l; j++)
      a(j, i) = a(i, j) = 0.0;
  }
}

// ###########################################################################
// reduction of a complex hermitian matrix to tridiagonal form
// Input:  a[n][n] is a std::real symmetric matrix
// Output: a[n][n] is orthogonal matrix Q effecting the transformation, it is
//                 an input for the subroutine tqli which finish diag.
//                 a_tridiag = Q^+ a Q
//         d[n] contains diagonal elements of the tridiagonal matrix,
//         e[n] contains LOWER subdiagonal elements, with e(1) = 0.
// based on tred2 from num.rec. and paper Mueller, num. math. 8, 72-92, 1966
// ###########################################################################

void tred2_c(MatrixXcd& a, int n, VectorXd& d, VectorXcd& e)
{
    int l, k, j, i;
    double scale, h;
    std::complex<double> g, f, hh;

    for (i = n - 1; i >= 1; i--)
    {
        l = i - 1;
        scale = 0.0;
        h = 0.0;
        if (l > 0)
        {
            for (k = 0; k <= l; k++)
                scale += abs(a(i, k));
            if (scale == 0.0)
                e(i) = a(i, l);
            else
            {
                for (k = 0; k <= l; k++)
                {
                    a(i, k) /= scale;
                    h += norm(a(i, k));
                }

                // my changes

                g = std::sqrt(h) * std::exp(std::complex<double>(0., arg(a(i, l))));
                h += std::sqrt(h) * abs(a(i, l));
                a(i, l) = a(i, l) + g;
                e(i) = -scale * g;

                f = std::complex<double>(0.0, 0.0);
                for (j = 0; j <= l; j++)
                {
                    a(j, i) = a(i, j) / h;
                    g = std::complex<double>(0.0, 0.0);
                    for (k = 0; k <= j; k++)
                        g += a(j, k) * conj(a(i, k));
                    for (k = j + 1; k <= l; k++)
                        g += conj(a(k, j)) * conj(a(i, k));
                    e(j) = g / h;
                    f += e(j) * a(i, j);
                }
                hh = f / (h + h);
                for (j = 0; j <= l; j++)
                {
                    f = conj(a(i, j));
                    e(j) = g = e(j) - hh * f;
                    for (k = 0; k <= j; k++)
                        a(j, k) -= (f * conj(e(k)) + g * a(i, k));
                }
            }
        }
        else
            e(i) = a(i, l);
        d(i) = h;
    }
    d(0) = 0.0;
    e(0) = std::complex<double>(0.0, 0.0);

    for (i = 0; i < n; i++)
    {
        l = i - 1;
        if (d(i))
        {
            for (j = 0; j <= l; j++)
            {
                g = std::complex<double>(0.0, 0.0);
                for (k = 0; k <= l; k++)
                    g += a(i, k) * a(k, j);
                for (k = 0; k <= l; k++)
                    a(k, j) -= g * conj(a(k, i));
            }
        }
        d(i) = std::real(a(i, i));
        a(i, i) = std::complex<double>(1.0, 0.0);
        for (j = 0; j <= l; j++)
            a(j, i) = a(i, j) = std::complex<double>(0.0, 0.0);
    }
}


// ###########################################################################
// suroutine for diagonalization of a std::real, symmetric, tridiagonal matrix.
// Input:  d[n] contains diagonal elements of the tridiagonal matrix,
//         e[n] contains subdiagonal elements with e(1) arbitrary.
//         z[n][n] should be identity matrix if eigenvectors of the tridiag
//                 matrix are desired,
//                 OR if matrix was previously reduced by tred2, then Z is the
//                 output matrix of that subroutin
// Output: d[n] contains eigenvalues,
//         e[n] is destroyed,
//         z[n][n] k-th column of z returns the normalized eigenvector
//                 corresponding to d[k]
// ###########################################################################

double SIGN(double r, double g)
{
  if (g < 0.0 ) return - fabs(r);
  else return fabs(r);
}

void tqli(VectorXd& d, VectorXd& e, int n, MatrixXd& z)
{
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

    for (i = 1; i < n; i++) e(i - 1) = e(i); // Adjusted indexing
    e(n - 1) = 0.0; // Adjusted indexing
    
    for (l = 0; l < n; l++) // Adjusted indexing
    {
        iter = 0;
        
        do
        {
            for (m = l; m < n - 1; m++) // Adjusted indexing
            {
                dd = fabs(d(m)) + fabs(d(m + 1)); // Adjusted indexing
                if ((float)(fabs(e(m)) + dd) == (float)dd) break;
            }

            if (m != l)
            {
                if (iter++ == 30) {
                    std::cout << "Too many iterations in tqli!" << std::endl;
                    std::cout.flush();
                    throw std::string("Too many iterations in tqli!");
                }
                g = (d(l + 1) - d(l)) / (2.0 * e(l)); // Adjusted indexing
                r = std::sqrt(g * g + 1.0);
                g = d(m) - d(l) + e(l) / (g + SIGN(r, g)); // Adjusted indexing
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--) // Adjusted indexing
                {
                    f = s * e(i); // Adjusted indexing
                    b = c * e(i); // Adjusted indexing
                    e(i + 1) = (r = std::sqrt(f * f + g * g)); // Adjusted indexing
                    if (r == 0.0)
                    {
                        d(i + 1) -= p; // Adjusted indexing
                        e(m) = 0.0; // Adjusted indexing
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d(i + 1) - p; // Adjusted indexing
                    r = (d(i) - g) * s + 2.0 * c * b; // Adjusted indexing
                    d(i + 1) = g + (p = s * r); // Adjusted indexing
                    g = c * r - b;

                    for (k = 0; k < n; k++) // Adjusted indexing
                    {
                        f = z(k, i + 1); // Adjusted indexing
                        z(k, i + 1) = s * z(k, i) + c * f; // Adjusted indexing
                        z(k, i) = c * z(k, i) - s * f; // Adjusted indexing
                    }

                }
                
                if (r == 0.0 && i >= l) continue;
                d(l) -= p; // Adjusted indexing
                e(l) = g; // Adjusted indexing
                e(m) = 0.0; // Adjusted indexing
            }

        } while (m != l);
    }
}



// ###########################################################################
// THE SAME AS ABOVE BUT z MATRIX IS COMPLEX DOUBLE
// suroutine for diagonalization of a std::real, symmetric, tridiagonal matrix.
// Input:  d[n] contains diagonal elements of the tridiagonal matrix,
//         e[n] contains subdiagonal elements with e(1) arbitrary.
//         z[n][n] should be identity matrix if eigenvectors of the tridiag
//                 matrix are desired,
//                 OR if matrix was previously reduced by tred2, then Z is the
//                 output matrix of that subroutin
// Output: d[n] contains eigenvalues,
//         e[n] is destroyed,
//         z[n][n] k-th column of z returns the normalized eigenvector
//                 corresponding to d[k]
// ###########################################################################

void tqli_c(VectorXd& d, VectorXd& e, int n, MatrixXcd& z)
{
    int m, l, iter, i, k;
    double s, r, p, g, f, dd, c, b;

    std::complex<double> fcom;

    for (i = 1; i < n; i++) e(i - 1) = e(i); // Changed indices to 0-based
    e(n - 1) = 0.0; // Changed index to 0-based

    for (l = 0; l < n; l++) // Changed loop start to 0-based
    {
        iter = 0;

        do
        {
            for (m = l; m < n - 1; m++) // Changed loop condition to 0-based
            {
                dd = fabs(d(m)) + fabs(d(m + 1)); // Changed indices to 0-based
                if ((float)(fabs(e(m)) + dd) == (float)dd) break;
            }

            if (m != l)
            {
                if (iter++ == 30) {
                    std::cout << "Too many iterations in tqli!" << std::endl;
                    std::cout.flush();
                    throw std::string("Too many iterations in tqli!");
                }
                g = (d(l + 1) - d(l)) / (2.0 * e(l)); // Changed indices to 0-based
                r = std::sqrt(g * g + 1.0);
                g = d(m) - d(l) + e(l) / (g + SIGN(r, g)); // Changed indices to 0-based
                s = c = 1.0;
                p = 0.0;

                for (i = m - 1; i >= l; i--) // Changed indices to 0-based
                {
                    f = s * e(i);
                    b = c * e(i);
                    e(i + 1) = (r = std::sqrt(f * f + g * g)); // Changed index to 0-based
                    if (r == 0.0)
                    {
                        d(i + 1) -= p; // Changed index to 0-based
                        e(m) = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d(i + 1) - p; // Changed index to 0-based
                    r = (d(i) - g) * s + 2.0 * c * b; // Changed index to 0-based
                    d(i + 1) = g + (p = s * r); // Changed index to 0-based
                    g = c * r - b;

                    for (k = 0; k < n; k++) // Changed loop start to 0-based
                    {
                        fcom = z(k, i + 1); // Changed index to 0-based
                        z(k, i + 1) = s * z(k, i) + c * fcom; // Changed index to 0-based
                        z(k, i) = c * z(k, i) - s * fcom; // Changed index to 0-based
                    }
                }

                if (r == 0.0 && i >= l) continue;
                d(l) -= p; // Changed index to 0-based
                e(l) = g; // Changed index to 0-based
                e(m) = 0.0;
            }

        } while (m != l);
    }
}


// ###########################################################################
// Diagonalization of a std::real symmetric matrix A
// OUTPUT:  eigenvalues in m[n], and eigenvectors in columns of A[n][n]
// diag (m) = Qdag A Q, where Q is A and original A is destroyed.
// ###########################################################################

void diag_r_s(MatrixXd &A, int n, VectorXd& m)
{
  
  VectorXd e(n);

// transforming A to tridiagonal form

  tred2(A, n, m, e);

// diagonalyzing the tridiag matrix

  tqli(m, e, n, A);

}


// ###########################################################################
// Diagonalization of a complex hermitian matrix A
// OUTPUT:  eigenvalues in m[n], and eigenvectors in columns of A[n][n]
// diag (m) = Qdag A Q, where Q is A and original A is destroyed.
// ###########################################################################

void diag_c_h(MatrixXcd& A, int n, VectorXd& m)
{

  VectorXd eR(n);
  VectorXcd e(n);

  double g;

// transforming A to tridiagonal form (still compl. herm.)

  tred2_c(A, n, m, e);

// transforming compl. herm. tridiag. to std::real symm. tridiag.

  g = 0.0;

  for(int j=0; j<n; j++)
  {
    g += arg(e(j));

    for(int i=0; i<n; i++)
    {
      A(i,j) = A(i,j) * exp(dcomplex(0.,g));
    }
    
    eR(j) = abs(e(j));
  }

  tqli_c(m, eR, n, A);

}




// ###########################################################################
// Bi-unitary transformation of a complex matrix A
// OUTPUT:  "eigenvalues" in m[n],
// diag (m) = VL^*  A  VR^T, where VL is A and original A is destroyed.
// this is a notation of T. Blazek, and Haber & Kane.
// ###########################################################################

void bi_un_tr(MatrixXcd& A, int n, VectorXd& m, MatrixXcd& VR)
{
    // temporary use UR to store Adag
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            VR(i, j) = conj(A(j, i));
        }
    }

    A = A * VR;

    // A is now A * Adag and so is hermitian
    diag_c_h(A, n, m);

    // eigenvalues of A Adag are m^2

    for (int i = 0; i < n; i++)
    {
        if (m(i) < 0.)
            m(i) = 0.;
        else
            m(i) = std::sqrt(m(i));
    }

    // ordering eigenvalues and eigenvectors
    // eigenvectors are in columns of A,
    // therefore exchanging eigenvalues = exchanging columns of A

    if (n == 3 && (m(0) > m(1) || m(1) > m(2) || m(0) > m(2)))
    {
        std::cout << "Ordering matrix!!!" << std::endl;
        order(n, m, A);
    }

    // calculating of Ydag * VL^T

    VR = VR * A;

    // VL^T is in A, after transposing VL is in A
    A = A.transpose().eval();

    // finishing calculating of VR, VR^T = Ydag * VL^T * m^{-1)
    for (int j = 0; j < n; j++)
    {
        if (m(j) == 0.)
        {
            for (int i = 0; i < n; i++)
            {
                VR(j, i) = std::complex<double>(0.0, 0.0);
                VR(i, j) = std::complex<double>(0.0, 0.0);
            }
            VR(j, j) = std::complex<double>(1.0, 0.0);
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                VR(i, j) = VR(i, j) / m(j);
            }
        }
    }

    VR = VR.transpose().eval();
}

double Diagonalize(double v, double tanb, double mu, VectorXd& g, VectorXcd& M,
        MatrixXcd& Yu, MatrixXcd& Yd, MatrixXcd& Ye,
        VectorXd& mYu, VectorXd& mYd, VectorXd& mYe, MatrixXcd& VCKM,
        MatrixXcd& TYu, MatrixXcd& TYd, MatrixXcd& TYe,
        MatrixXcd& MQ2, MatrixXcd& ML2, MatrixXcd& Mu2,
        MatrixXcd& Md2, MatrixXcd& Me2,
        MatrixXcd& mQu, MatrixXcd& mQd,
        VectorXd& mMu, VectorXd& mMd, VectorXd& mMe,
        MatrixXcd& GuL, MatrixXcd& GdL, MatrixXcd& GeL,
        MatrixXcd& GuR, MatrixXcd& GdR, MatrixXcd& GeR,
        MatrixXcd& GvL,
        VectorXd& mMch, MatrixXcd& U, MatrixXcd& V,
        VectorXd& mMnt, MatrixXcd&Nt, MatrixXd& Nt_mix,
        double MZ, double MW, double sw2, double mA, double mHp, double mH,
        double mh, double c2a, double s2a,
        double& mstop_1, double& mstop_2, double& s2stop, double& c2stop,
        double& msbot_1, double& msbot_2, double& s2sbot, double& c2sbot,
        double& mstau_1, double& mstau_2, double& s2stau, double& c2stau, double& Mg)
{
    
    double MEW = 91.188;
    const double pi = 3.141592653589793;
    double P;
    P = 1.;
    double vu, vd;

    MatrixXcd VLu(3,3), VLd(3,3), VLe(3,3), VLn(3,3);
    MatrixXcd VRu(3,3), VRd(3,3), VRe(3,3), VRn(3,3);


    // y - hypercharge
    double yq = 1./3., yu = -4./3., yd = 2./3., yl = -1., ye = 2.;

    double DuLL, DuRR, DdLL, DdRR, DeLL, DeRR, DvLL;

    // sparticle mass matrices
    MatrixXcd M2u(6,6), M2d(6,6), M2e(6,6), M2v(3,3), pom(3,3);

    MatrixXd Ntr(4,4);
    VectorXcd PNt(4);

    for(int i=0; i<4; i++) PNt(i) = dcomplex(1.,0);

    MatrixXcd mMnpom(3,3), PVLn(3,3);


    //Diagonalization of Yu, ...;
    VLu = Yu;
    VLd = Yd;
    VLe = Ye;

    bi_un_tr(VLu, 3, mYu, VRu);
    bi_un_tr(VLd, 3, mYd, VRd);
    bi_un_tr(VLe, 3, mYe, VRe);

    // CKM matrix
    VCKM = (VLu * VLd.adjoint()).eval();

    vd = v / std::sqrt(2. + 2.* tanb * tanb);
    vu = v / std::sqrt(2. + 2. /(tanb*tanb));

    //cout<<"End of Diagonalization of Yu, ..."<<"\n";
    //cout<<"v = "<<v<<"  vu = "<<vu<<"  vd = "<<vd<<"\n";
    /*
    cout<<"just diagonal Yukawas"<<endl;
    cout<<"mYu: "<<mYu<<"\n";
    cout<<"mYd: "<<mYd<<"\n";
    cout<<"mYe: "<<mYe<<"\n";
    */
    // fermion masses
    for (int i = 0; i < 3; i++)
    {
      mYu(i) = mYu(i)*vu;
      mYd(i) = mYd(i)*vd;
      mYe(i) = mYe(i)*vd;
    }

  // making MQ2 matrices hermitian by throwing away diagonal phases

    for (int i = 0; i <= 2; i++)
    {
      MQ2(i, i) = dcomplex(std::real(MQ2(i, i)), 0.);
      ML2(i, i) = dcomplex(std::real(ML2(i, i)), 0.);
      Mu2(i, i) = dcomplex(std::real(Mu2(i, i)), 0.);
      Md2(i, i) = dcomplex(std::real(Md2(i, i)), 0.);
      Me2(i, i) = dcomplex(std::real(Me2(i, i)), 0.);
    }

  //  cout<<"Generating squark and slepton matrices"<<"\n";

    DuLL = (vu*vu - vd*vd)*(-g[1]*g[1] + yq*g[0]*g[0]*3./5.)/4.;
    DuRR = (vu*vu - vd*vd)* yu * g[0]*g[0]*(3./5.)/4.;
    DdLL = (vu*vu - vd*vd)*( g[1]*g[1] + yq*g[0]*g[0]*3./5.)/4.;
    DdRR = (vu*vu - vd*vd)* yd * g[0]*g[0]*(3./5.)/4.;
    DeLL = (vu*vu - vd*vd)*( g[1]*g[1] + yl*g[0]*g[0]*3./5.)/4.;
    DeRR = (vu*vu - vd*vd)* ye * g[0]*g[0]*(3./5.)/4.;
    DvLL = (vu*vu - vd*vd)*(-g[1]*g[1] + yl*g[0]*g[0]*3./5.)/4.;

    M2u.block<3, 3>(0, 0) = MQ2 + (vu * vu) * Yu.conjugate() * Yu.transpose().eval();
    M2u.block<3, 3>(0, 3) = (vu * TYu.conjugate() - (mu * vd) * Yu.conjugate());
    M2u.block<3, 3>(3, 0) = (vu * TYu - (mu * vd) * Yu);
    M2u.block<3, 3>(3, 3) = Mu2 + (vu * vu) * Yu.transpose().eval() * Yu.conjugate();

    M2d.block<3, 3>(0, 0) = MQ2 + (vd * vd) * Yd.conjugate() * Yd.transpose().eval();
    M2d.block<3, 3>(0, 3) = (vd * TYu.conjugate() - (mu * vu) * Yd.conjugate());
    M2d.block<3, 3>(3, 0) = (vd * TYd - (mu * vu) * Yd);
    M2d.block<3, 3>(3, 3) = Md2 + (vd * vd) * Yd.transpose().eval() * Yd.conjugate();

    M2e.block<3, 3>(0, 0) = ML2 + (vd * vd) * Ye.conjugate() * Ye.transpose().eval();
    M2e.block<3, 3>(0, 3) = (vd * TYe.conjugate() - (mu * vu) * Ye.conjugate());
    M2e.block<3, 3>(3, 0) = (vd * TYe - (mu * vu) * Ye);
    M2e.block<3, 3>(3, 3) = Me2 + (vd * vd) * Ye.transpose().eval() * Ye.conjugate();


    for (int i = 0; i < 3; i++) {
        M2u(i, i) += DuLL;
        M2u(i + 3, i + 3) += DuRR;

        M2d(i, i) += DdLL;
        M2d(i + 3, i + 3) += DdRR;

        M2e(i, i) += DeLL;
        M2e(i + 3, i + 3) += DeRR;
    }


    // ROTATE TO CKM AND PMNS BASIS !!!
    TYu = ((VLu).conjugate() * TYu * (VRu).transpose()).eval();
    TYd = ((VLd).conjugate() * TYd * (VRd).transpose()).eval();
    TYe = ((VLe).conjugate() * TYe * (VRe).transpose()).eval();

    MQ2 = (VLd * MQ2 * (VLd).adjoint()).eval();             // V_d^+ MQ2 V_d   in SLHA2 language (Eq.(13) in arXiv:0801.0045)
    Mu2 = (VRu * (Mu2).transpose() * (VRu).adjoint()).eval();  // U_u^+ Mu2^T U_u in SLHA2 language (Eq.(13) in arXiv:0801.0045)
    Md2 = (VRd * (Md2).transpose() * (VRd).adjoint()).eval();  // U_d^+ MQ2^T U_d in SLHA2 language (Eq.(13) in arXiv:0801.0045)

    ML2 = (VLe * ML2 * (VLe).adjoint()).eval();             // V_d^+ MQ2 V_d   in SLHA2 language (Eq.(13) in arXiv:0801.0045)
    Me2 = (VRe * (Me2).transpose() * (VRe).adjoint()).eval();  // U_d^+ MQ2^T U_d in SLHA2 language (Eq.(13) in arXiv:0801.0045)
    // Analogous term for rh sneutrino missing!


  /*
    cout<<"M2u matrix: "<<M2u<<"\n";
    cout<<"M2d matrix: "<<M2d<<"\n";
    cout<<"M2e matrix: "<<M2e<<"\n";
    cout<<"M2v matrix: "<<M2v<<"\n";
  */

    diag_c_h(M2u, 6, mMu);
    diag_c_h(M2d, 6, mMd);
    diag_c_h(M2e, 6, mMe);


  // the eigenvectors are in M2u, ...

    M2u = ((M2u).adjoint()).eval();
    M2d = ((M2d).adjoint()).eval();
    M2e = ((M2e).adjoint()).eval();
    M2v = ((M2v).adjoint()).eval();

  //  cout<<"M2u = "<<mMu<<endl;

    for (int i = 0; i <= 5; i++)
    {
      for (int j = 0; j <= 2; j++)
      {
        GuL(i, j) = M2u(i, 0)*conj(VLu(j, 0))+M2u(i, 1)*conj(VLu(j, 1))+M2u(i, 2)*conj(VLu(j, 2));
        GdL(i, j) = M2d(i, 0)*conj(VLd(j, 0))+M2d(i, 1)*conj(VLd(j, 1))+M2d(i, 2)*conj(VLd(j, 2));
        GeL(i, j) = M2e(i, 0)*conj(VLe(j, 0))+M2e(i, 1)*conj(VLe(j, 1))+M2e(i, 2)*conj(VLe(j, 2));

        GuR(i, j) = M2u(i, 3)*conj(VRu(j, 0))+M2u(i, 4)*conj(VRu(j, 1))+M2u(i, 5)*conj(VRu(j, 2));
        GdR(i, j) = M2d(i, 3)*conj(VRd(j, 0))+M2d(i, 4)*conj(VRd(j, 1))+M2d(i, 5)*conj(VRd(j, 2));
        GeR(i, j) = M2e(i, 3)*conj(VRe(j, 0))+M2e(i, 4)*conj(VRe(j, 1))+M2e(i, 5)*conj(VRe(j, 2));
      }
    }

  // from diag. we have masses squared

    for (int i = 0; i <= 2; i++)
    {
      if (mMu[i]>=0.) mMu[i] = std::sqrt(mMu[i]);
      else P = P*10.;

      if (mMd[i]>=0.) mMd[i] = std::sqrt(mMd[i]);
      else P = P*10.;

      if (mMe[i]>=0.) mMe[i] = std::sqrt(mMe[i]);
      else P = P*10.;

      if (mMu[i+3]>=0.) mMu[i+3] = std::sqrt(mMu[i+3]);
      else P = P*10.;

      if (mMd[i+3]>=0.) mMd[i+3] = std::sqrt(mMd[i+3]);
      else P = P*10.;

      if (mMe[i+3]>=0.) mMe[i+3] = std::sqrt(mMe[i+3]);
      else P = P*10.;
    }

    if (P!=1.) return P;

  //  cout<<"End of Diagonalization of squarks, ..."<<"\n";


  // trilinear matrices replaced by new matrices TYu = VLu^* TYu VRu^T
  /*  yu = (VLu).conjugate() * yu * (VRu).transpose().eval();
    yd = (VLd).conjugate() * yd * (VRd).transpose().eval();
    ye = (VLe).conjugate() * ye * (VRe).transpose().eval();
  */


  /*
    cout<<"TYu: "<<TYu<<"\n";
    cout<<"TYd: "<<TYd<<"\n";
    cout<<"TYe: "<<TYe<<"\n";
  */
  // no generation mixing, just 3rd gen. L-R mixing approximation for ewsb !!

  // Dterms are not universal for up and down!
  // for sleptons ML2 modified for selectrons, for sneutrino not needed!!!
  // add D-terms to 1st two gens.


  //  mQu = MQ2;
  //  mQd = MQ2;

    for (int i = 0; i < 2; i++)
    {
      mQu(i,i) = MQ2(i,i) + DuLL;
      mQd(i,i) = MQ2(i,i) + DdLL;
      Mu2(i,i) = Mu2(i,i) + DuRR;
      Md2(i,i) = Md2(i,i) + DdRR;
      ML2(i,i) = ML2(i,i) + DeLL;
      Me2(i,i) = Me2(i,i) + DeRR;
    }
  /*
    cout<<"mQu: "<<mQu<<"\n";
    cout<<"mQd: "<<mQd<<"\n";
    cout<<"Mu2: "<<Mu2<<"\n";
    cout<<"Md2: "<<Md2<<"\n";
  */
    double a,b,d;

    double sstop2, cstop2, sstop, cstop;

    MatrixXd G_stop(2,2), M_stop(2,2);

    a = std::real(MQ2(2,2)) + mYu(2)*mYu(2) + DuLL;
    b = vu*std::real(TYu(2,2)) - mu*mYu(2)/tanb;
    d = std::real(Mu2(2,2)) + mYu(2)*mYu(2) + DuRR;
  /*
    cout<<"a = "<<a<<endl;
    cout<<"b = "<<b<<endl;
    cout<<"d = "<<d<<endl;
  */

    mstop_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
    mstop_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b) ;

  //  cout<<"mstop1 ="<<mstop_1<<endl;
  //  cout<<"mstop2 = "<<mstop_2<<endl;

    if (mstop_2 >=0.) mstop_2 = std::sqrt(mstop_2);
    else P = P*10.;

    if (mstop_1 >=0.) mstop_1 = std::sqrt(mstop_1);
    else P = P*10.;

    if (P!=1.) return P;

    s2stop = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
    c2stop = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);

  //  double msbot_1, msbot_2, s2sbot, c2sbot, ssbot2, csbot2;

    a = std::real(MQ2(2,2)) + mYd(2)*mYd(2) + DdLL;
    b = vd*std::real(TYd(2,2)) - mu*mYd(2)*tanb;
    d = std::real(Md2(2,2)) + mYd(2)*mYd(2) + DdRR;

    msbot_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
    msbot_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b) ;

    if (msbot_2 >=0.) msbot_2 = std::sqrt(msbot_2);
    else P = P*10.;

    if (msbot_1 >=0.) msbot_1 = std::sqrt(msbot_1);
    else P = P*10.;

    if (P!=1.) return P;

    s2sbot = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
    c2sbot = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);

  //  double mstau_1, mstau_2, s2stau, c2stau, sstau2, cstau2;

    a = std::real(ML2(2,2)) + mYe(2)*mYe(2) + DeLL;
    b = vd*std::real(TYe(2,2)) - mu*mYe(2)*tanb;
    d = std::real(Me2(2,2)) + mYe(2)*mYe(2) + DeRR;

    mstau_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
    mstau_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b) ;

    if (mstau_2 >=0.) mstau_2 = std::sqrt(mstau_2);
    else P = P*10.;

    if (mstau_1 >=0.) mstau_1 = std::sqrt(mstau_1);
    else P = P*10.;

    if (P!=1.) return P;

    s2stau = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
    c2stau = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);

  // chargino mass matrix
    U(0, 0) = dcomplex(std::real(M(1)), 0.);
    U(0, 1) = g(1) * vu;
    U(1, 0) = g(1) * vd;
    U(1, 1) = mu;

  //  cout<<"U =" << U;

    MatrixXcd Ch_mass(2,2), Ch_diag(2,2);
    Ch_mass = U;

    bi_un_tr(U, 2, mMch, V);

    MatrixXd U_mix(2,2), V_mix(2,2);
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 2; j++)
        {
        U_mix(i,j) = std::real(U(i,j));
        V_mix(i,j) = std::real(V(j,i));
        }
     V = (V).conjugate();
  /*
    cout<<"mMch = "<<mMch<<endl;
    cout<<"U_mix = "<<U_mix<<endl;
    cout<<"V_mix = "<<V_mix<<endl;
  */
  //This part is to check if the matrix U and V diagonalize the chargino mass matrix.
  /*
   for(int i=0; i< 2; i++)
   {
     for(int j=0; j<2; j++)
     {
      Ch_diag(i, j) = 0;
      for(int k=0; k<2; k++)
      {
       for(int l=0; l<2; l++)
       {
        Ch_diag(i, j) = Ch_diag(i, j) + (U_mix(i, k)*Ch_mass(k, l)*V_mix(l, j));
       }
      }
     }
   }
   cout<<"Ch_diag = "<<Ch_diag;
  */
  // neutralino mass matrix
    Ntr(0, 0) = std::real(M[0]);
    Ntr(0, 1) = 0.;
    Ntr(0, 2) = - g[0]*std::sqrt(3./5) * vd/std::sqrt(2.);
    Ntr(0, 3) = g[0]*std::sqrt(3./5) * vu/std::sqrt(2.);

    Ntr(1, 1) = std::real(M[1]);
    Ntr(1, 2) = g[1] * vd/std::sqrt(2.);
    Ntr(1, 3) = - g[1] * vu/std::sqrt(2.);

    Ntr(2, 2) = 0.;
    Ntr(2, 3) = - mu;

    Ntr(3, 3) = 0.;

    for (int i = 1; i <= 3; i++)
    {
      for (int j = 0; j < i; j++)
      {
        Ntr(i, j) = Ntr(j, i);
      }
    }

  // cout<<"Ntr = "<<Ntr<<endl;

  MatrixXcd Ntr_zal(4,4);
    MatrixXd Ntr_mass(4,4);
    for (int i = 0; i < 4; i++)
      for (int j =0; j < 4; j++)
        {
        Ntr_zal(i,j) = Ntr(i,j);
        Ntr_mass(i,j) = Ntr(i,j);
        }

    diag_r_s(Ntr, 4, mMnt);
  //  cout<<"Ntr = " <<Ntr<<endl;
  // result is Ntr^T
    Nt_mix = Ntr;
    Ntr = (Ntr).transpose().eval();
  //  cout<<"After transposing"<<endl;
  //  cout<<"Ntr = "<<Ntr;

  //This part is to check if the matrix Nt_mix and Ntr diagonalize the mass matrix.
  /*
    MatrixXd Nt_diag(4,4);
   for(int i=0; i< 4; i++)
   {
     for(int j=0; j<4; j++)
     {
      Nt_diag(i, j) = 0;
      for(int k=0; k<4; k++)
      {
       for(int l=0; l<4; l++)
       {
        Nt_diag(i, j) = Nt_diag(i, j) + (Ntr(i, k)*Ntr_mass(k, l)*Nt_mix(l, j));
       }
      }
     }
   }
   cout<<"Nt_diag = "<<Nt_diag;
  */
  //  Nt_mix = (Nt_mix).transpose().eval(); //This transposed quantity is what you want in the SLHA output.
    for (int i = 0; i < 4; i++)
    {
      if (mMnt(i) < 0.)
      {
        mMnt(i) = -mMnt(i);
  // forming P^*
        PNt(i) = dcomplex(0.,-1.);
      }
    }


    for (int i = 0; i < 4; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        Nt(i,j) = PNt(i)*Ntr(i,j);
      }
    }

  //just checking
  /*
  MatrixXcd Njunk(4,4);

    for (int i = 1; i <= 4; i++)
      for (int j = 1; j <= 4; j++)
      {
        Njunk(i,j) = dcomplex(0.,0.);
        for (int k = 1; k <= 4; k++)
          for (int l = 1; l <= 4; l++)
            Njunk(i,j) = Njunk(i,j) + conj(Nt(i,k))*Ntr_zal(k,l)*conj(Nt(j,l));
      }

  //cout<<"mMnt = "<<mMnt<<"\n";
  //cout<<"Nt = "<<Nt<<"\n";
  //cout<<"mMnt = "<<Njunk<<"\n";
  */
  /*
  cout<<"mMnt = "<<mMnt<<"\n";
  cout<<"Nt = "<<Nt<<"\n";
  */
  // Gluino pole mass!!!!!!!!!!

  // another calculation
  /*
    Mg = (15. + 18.*log(MEW/std::real(M(3))))*g(3)*g(3)/(16.*pi*pi);

   for (int a = 1; a <= 6; a++)
      Mg = Mg - ( B1(std::real(M(3)), 0., mMu(a)) +
           B1(std::real(M(3)), 0., mMd(a)) )*g(3)*g(3)/(16.*pi*pi);

    Mg = std::real(M(3))/(1. - Mg);

    cout<<"M(3) = "<<std::real(M(3))<<endl;
    cout<<"Mg = "<<Mg<<endl;
  */
  // start from beginning!!!!!!!!!!!!!!!!!!!!!

    Mg = std::real(M(2)) - ( - std::real(M(2))
         * (15. + 18.*std::log(MEW/std::real(M(2)))) )*g(2)*g(2)/(16.*pi*pi);

    for (int a = 0; a < 6; a++)
      Mg = Mg + std::real(M(2)) * ( B1(std::real(M(2)), 0., mMu(a)) +
           B1(std::real(M(2)), 0., mMd(a)) )*g(2)*g(2)/(16.*pi*pi);

  //  cout<<"Mg = "<<Mg<<endl;


  //cout<<"Yu = "<<Yu<<endl;

  //cout<<"thresh_flag = "<<thresh_flag<<endl;
  // threshold corrections to fermion masses
  ThresholdCorrectionsFermions(v, tanb, g, M, mYu, mYd, mYe, VCKM, Yu, Yd, Ye,
     mMu, mMd, mMe, GuL, GdL, GeL, GuR, GdR, GeR, GvL, mMch, U, V, mMnt, Nt,
     MZ, MW, sw2, mA, mHp, mH, mh, c2a, s2a);
      // after threshold corrections Yu, ... contain new Yukawa matrices (multiplied by vevs)
      // in the basis where original Yukawa matrices were diagonal!!!

      // fermion yukawas after corrections

      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          Yu(i,j) = Yu(i,j)/vu;
          Yd(i,j) = Yd(i,j)/vd;
          Ye(i,j) = Ye(i,j)/vd;
        }
      }
      /*
      cout<<"After corrections: "<<"\n";
      cout<<"Yu matrix: "<<Yu<<"\n";
      cout<<"Yd matrix: "<<Yd<<"\n";
      cout<<"Ye matrix: "<<Ye<<"\n";


      for (int i = 1; i <= 3; i++)
      {
      for (int j = 1; j <= 3; j++)
      {
        Yu(i,j) = Yu(i,j)*vu;
        Yd(i,j) = Yd(i,j)*vd;
        Ye(i,j) = Ye(i,j)*vd;
      }
      }
      */
      //cout<<"RE-Diagonalization of Yu, ..."<<"\n";

      VLu = Yu;
      VLd = Yd;
      VLe = Ye;

      bi_un_tr(VLu, 3, mYu, VRu);
      bi_un_tr(VLd, 3, mYd, VRd);
      bi_un_tr(VLe, 3, mYe, VRe);

      // fermion masses
      for (int i = 0; i < 3; i++) {
        mYu(i) = mYu(i)*vu;
        mYd(i) = mYd(i)*vd;
        mYe(i) = mYe(i)*vd;
      }

      // new CKM matrix
      /*
      cout<<"VLd = "<<VLd<<endl;
      cout<<" (VLd).adjoint() = "<< (VLd).adjoint()<<endl;
      cout<<"VLd = "<<VLd<<endl;
      cout<<"VCKM = "<<VLu * VCKM * (VLd).adjoint()<<endl;
      */

      VCKM = VLu * VCKM * (VLd).adjoint();

      //  cout<<"VCKMcor = "<<VCKM<<endl;
      /*
      for (int i = 1; i <= 3; i++)
      {
      for (int j = 1; j <= 3; j++)
      {
        ABSVCKM(i,j) = abs(VCKM(i,j));
      }
      }
      */
      //  cout<<"ABSVCKMcor matrix:"<<ABSVCKM<<"\n";

      /*
      cout<<"Yu matrix: "<<Yu<<"\n";
      cout<<"Yd matrix: "<<Yd<<"\n";
      cout<<"Ye matrix: "<<Ye<<"\n";


      cout<<"mYu: "<<mYu<<"\n";
      cout<<"mYd: "<<mYd<<"\n";
      cout<<"mYe: "<<mYe<<"\n";



      // fermion yukawas after corrections
      for (int i = 1; i <= 3; i++)
      {
        mYu(i) = mYu(i)/vu;
        mYd(i) = mYd(i)/vd;
        mYe(i) = mYe(i)/vd;
      }

      cout<<"diagonal Yukawas after corrections"<<endl;
      cout<<"mYu: "<<mYu<<"\n";
      cout<<"mYd: "<<mYd<<"\n";
      cout<<"mYe: "<<mYe<<"\n";

      // fermion yukawas after corrections
      for (int i = 1; i <= 3; i++)
      {
        mYu(i) = mYu(i)*vu;
        mYd(i) = mYd(i)*vd;
        mYe(i) = mYe(i)*vd;
      }

     cout<<"diagonal Yukawas after corrections"<<endl;
      cout<<"mYu: "<<mYu<<"\n";
      cout<<"mYd: "<<mYd<<"\n";
      cout<<"mYe: "<<mYe<<"\n";
      */

      //Repeating the susy spectrum calculation after threshold corrections:
      //  cout<<"Generating squark and slepton matrices"<<"\n";

      DuLL = (vu*vu - vd*vd)*(-g[1]*g[1] + yq*g[0]*g[0]*3./5.)/4.;
      DuRR = (vu*vu - vd*vd)* yu * g[0]*g[0]*(3./5.)/4.;
      DdLL = (vu*vu - vd*vd)*( g[1]*g[1] + yq*g[0]*g[0]*3./5.)/4.;
      DdRR = (vu*vu - vd*vd)* yd * g[0]*g[0]*(3./5.)/4.;
      DeLL = (vu*vu - vd*vd)*( g[1]*g[1] + yl*g[0]*g[0]*3./5.)/4.;
      DeRR = (vu*vu - vd*vd)* ye * g[0]*g[0]*(3./5.)/4.;
      DvLL = (vu*vu - vd*vd)*(-g[1]*g[1] + yl*g[0]*g[0]*3./5.)/4.;

    
    M2u.block<3, 3>(0, 0) = (VCKM * MQ2 * VCKM.adjoint() + (vu * vu) * Yu.conjugate() * Yu.transpose()).eval();
    M2u.block<3, 3>(0, 3) = (vu * TYu.conjugate() - (mu * vd) * Yu.conjugate()).eval();
    M2u.block<3, 3>(3, 0) = vu * TYu - (mu * vd) * Yu;
    M2u.block<3, 3>(3, 3) = (Mu2 + (vu * vu) * Yu.transpose().eval() * Yu.conjugate()).eval();

    M2d.block<3, 3>(0, 0) = (MQ2 + (vd * vd) * Yd.conjugate() * Yd.transpose()).eval();
    M2d.block<3, 3>(0, 3) = (vd * TYu.conjugate() - (mu * vu) * Yd.conjugate()).eval();
    M2d.block<3, 3>(3, 0) = vd * TYd - (mu * vu) * Yd;
    M2d.block<3, 3>(3, 3) = (Md2 + (vd * vd) * Yd.transpose().eval() * Yd.conjugate()).eval();

    M2e.block<3, 3>(0, 0) = (ML2 + (vd * vd) * Ye.conjugate() * Ye.transpose()).eval();
    M2e.block<3, 3>(0, 3) = (vd * TYe.conjugate() - (mu * vu) * Ye.conjugate()).eval();
    M2e.block<3, 3>(3, 0) = vd * TYe - (mu * vu) * Ye;
    M2e.block<3, 3>(3, 3) = (Me2 + (vd * vd) * Ye.transpose().eval() * Ye.conjugate()).eval();

    for (int i = 0; i < 3; i++) {
        M2u(i, i) += DuLL;
        M2u(i + 3, i + 3) += DuRR;

        M2d(i, i) += DdLL;
        M2d(i + 3, i + 3) += DdRR;

        M2e(i, i) += DeLL;
        M2e(i + 3, i + 3) += DeRR;
    }


      // ROTATE TO CKM AND PMNS BASIS !!!
      TYu = ((VLu).conjugate() * TYu * (VRu).transpose()).eval();
      TYd = ((VLd).conjugate() * TYd * (VRd).transpose()).eval();
      TYe = ((VLe).conjugate() * TYe * (VRe).transpose()).eval();

      MQ2 = (VLd * MQ2 * (VLd).adjoint()).eval();             // V_d^+ MQ2 V_d   in SLHA2 language (Eq.(13) in arXiv:0801.0045)
      Mu2 = (VRu * (Mu2).transpose().eval() * (VRu).adjoint()).eval();  // U_u^+ Mu2^T U_u in SLHA2 language (Eq.(13) in arXiv:0801.0045)
      Md2 = (VRd * (Md2).transpose().eval() * (VRd).adjoint()).eval();  // U_d^+ MQ2^T U_d in SLHA2 language (Eq.(13) in arXiv:0801.0045)

      ML2 = (VLe * ML2 * (VLe).adjoint()).eval();             // V_d^+ MQ2 V_d   in SLHA2 language (Eq.(13) in arXiv:0801.0045)
      Me2 = (VRe * (Me2).transpose().eval() * (VRe).adjoint()).eval();  // U_d^+ MQ2^T U_d in SLHA2 language (Eq.(13) in arXiv:0801.0045)
      // Analogous term for rh sneutrino missing!

      /*
      cout<<"M2u matrix: "<<M2u<<"\n";
      cout<<"M2d matrix: "<<M2d<<"\n";
      cout<<"M2e matrix: "<<M2e<<"\n";
      cout<<"M2v matrix: "<<M2v<<"\n";
      */

      diag_c_h(M2u, 6, mMu);
      diag_c_h(M2d, 6, mMd);
      diag_c_h(M2e, 6, mMe);

      // the eigenvectors are in M2u, ...

      M2u = (M2u).adjoint().eval();
      M2d = (M2d).adjoint().eval();
      M2e = (M2e).adjoint().eval();

      //  cout<<"M2u = "<<mMu<<endl;

      for (int i = 0; i <= 5; i++) {
        for (int j = 0; j <= 2; j++) {
          GuL(i, j) = M2u(i, 0)*conj(VLu(j, 0))+M2u(i, 1)*conj(VLu(j, 1))+M2u(i, 2)*conj(VLu(j, 2));
          GdL(i, j) = M2d(i, 0)*conj(VLd(j, 0))+M2d(i, 1)*conj(VLd(j, 1))+M2d(i, 2)*conj(VLd(j, 2));
          GeL(i, j) = M2e(i, 0)*conj(VLe(j, 0))+M2e(i, 1)*conj(VLe(j, 1))+M2e(i, 2)*conj(VLe(j, 2));

          GuR(i, j) = M2u(i, 3)*conj(VRu(j, 0))+M2u(i, 4)*(VRu(j, 1))+M2u(i, 5)*conj(VRu(j, 2));
          GdR(i, j) = M2d(i, 3)*conj(VRd(j, 0))+M2d(i, 4)*(VRd(j, 1))+M2d(i, 5)*conj(VRd(j, 2));
          GeR(i, j) = M2e(i, 3)*conj(VRe(j, 0))+M2e(i, 4)*conj(VRe(j, 1))+M2e(i, 5)*conj(VRe(j, 2));
        }
      }

      // from diag. we have masses squared

      for (int i = 0; i <= 2; i++) {
        if (mMu[i]>=0.) mMu[i] = std::sqrt(mMu[i]);
        else P = P*10.;

        if (mMd[i]>=0.) mMd[i] = std::sqrt(mMd[i]);
        else P = P*10.;

        if (mMe[i]>=0.) mMe[i] = std::sqrt(mMe[i]);
        else P = P*10.;

        if (mMu[i+3]>=0.) mMu[i+3] = std::sqrt(mMu[i+3]);
        else P = P*10.;

        if (mMd[i+3]>=0.) mMd[i+3] = std::sqrt(mMd[i+3]);
        else P = P*10.;

        if (mMe[i+3]>=0.) mMe[i+3] = std::sqrt(mMe[i+3]);
        else P = P*10.;
      }

      if (P!=1.) return P;

      //  cout<<"End of Diagonalization of squarks, ..."<<"\n";

      /*
      cout<<"GuL matrix:"<<GuL<<"\n";
      cout<<"GuR matrix:"<<GuR<<"\n";
      cout<<"GdL matrix:"<<GdL<<"\n";
      cout<<"GdR matrix:"<<GdR<<"\n";
      cout<<"GeL matrix:"<<GeL<<"\n";
      cout<<"GeR matrix:"<<GeR<<"\n";
      cout<<"GvL matrix:"<<GvL<<"\n";

      cout<<"mMu: "<<mMu<<"\n";
      cout<<"mMd: "<<mMd<<"\n";
      cout<<"mMe: "<<mMe<<"\n";
      cout<<"mMv: "<<mMv<<"\n";
      */

      // trilinear matrices replaced by new matrices TYu = VLu^* TYu VRu^T
      /*
      yu = (VLu).conjugate() * yu * (VRu).transpose().eval();
      yd = (VLd).conjugate() * yd * (VRd).transpose().eval();
      ye = (VLe).conjugate() * ye * (VRe).transpose().eval();
      */


      /*
      cout<<"TYu: "<<TYu<<"\n";
      cout<<"TYd: "<<TYd<<"\n";
      cout<<"TYe: "<<TYe<<"\n";
      */

      // no generation mixing, just 3rd gen. L-R mixing approximation for ewsb !!

      // Dterms are not universal for up and down!
      // for sleptons ML2 modified for selectrons, for sneutrino not needed!!!
      // add D-terms to 1st two gens.


      //  mQu = MQ2;
      //  mQd = MQ2;

      for (int i = 0; i < 2; i++) {
        mQu(i,i) = MQ2(i,i) + DuLL;
        mQd(i,i) = MQ2(i,i) + DdLL;
        Mu2(i,i) = Mu2(i,i) + DuRR;
        Md2(i,i) = Md2(i,i) + DdRR;
        ML2(i,i) = ML2(i,i) + DeLL;
        Me2(i,i) = Me2(i,i) + DeRR;
      }
      /*
      cout<<"mQu: "<<mQu<<"\n";
      cout<<"mQd: "<<mQd<<"\n";
      cout<<"Mu2: "<<Mu2<<"\n";
      cout<<"Md2: "<<Md2<<"\n";
      */

      a = std::real(MQ2(2,2)) + mYu(2)*mYu(2) + DuLL;
      b = vu*std::real(TYu(2,2)) - mu*mYu(2)/tanb;
      d = std::real(Mu2(2,2)) + mYu(2)*mYu(2) + DuRR;
      /*
      cout<<"a = "<<a<<endl;
      cout<<"b = "<<b<<endl;
      cout<<"d = "<<d<<endl;
      */
      mstop_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
      mstop_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b) ;

      //cout<<"mstop1 ="<<mstop_1<<endl;
      //cout<<"mstop2 = "<<mstop_2<<endl;

      if (mstop_2 >=0.) mstop_2 = std::sqrt(mstop_2);
      else P = P*10.;

      if (mstop_1 >=0.) mstop_1 = std::sqrt(mstop_1);
      else P = P*10.;

      if (P!=1.) return P;

      s2stop = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
      c2stop = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);

      /*
      sstop2 = 0.5 - 0.5*c2stop;
      cstop2 = 1. - sstop2;

      sstop = std::sqrt(sstop2)*sstop2/fabs(sstop2);
      cstop = std::sqrt(cstop2);

      G_stop(1,1) = cstop;
      G_stop(1,2) = sstop;
      G_stop(2,1) = -sstop;
      G_stop(2,2) = cstop;

      M_stop(1,1) = a;
      M_stop(1,2) = b;
      M_stop(2,1) = b;
      M_stop(2,2) = d;

      //Check if M_stop is diagonal
      cout<<"G*M*tr(G) = "<<G_stop*M_stop*(G_stop).transpose().eval()<<"\n";

      cout<<"mt_1^2 = "<<mstop_1*mstop_1<<"\n";
      cout<<"mt_2^2 = "<<mstop_2*mstop_2<<"\n";

      cout<<"std::sqrt(a*ctop2 + d*stop2 + b*s2top) = "<<std::sqrt(a*cstop2 + d*sstop2 + b*s2stop)<<"\n";
      cout<<"std::sqrt(a*stop2 + d*ctop2 - b*s2top) = "<<std::sqrt(a*sstop2 + d*cstop2 - b*s2stop)<<"\n";
      cout<<"b*c2top - 0.5*(a-d)*s2top = "<<b*c2stop - 0.5*(a-d)*s2stop<<"\n";
      cout<<"ctop2 = "<<cstop2<<"\n";
      cout<<"std::sqrt(ctop2) = "<<std::sqrt(cstop2)<<"\n";
      cout<<"stop2 = "<<sstop2<<"\n";
      cout<<"std::sqrt(stop2) = "<<std::sqrt(sstop2)<<"\n";
      */

      //double msbot_1, msbot_2, s2sbot, c2sbot,
      //double ssbot2, csbot2;

      a = std::real(MQ2(2,2)) + mYd(2)*mYd(2) + DdLL;
      b = vd*std::real(TYd(2,2)) - mu*mYd(2)*tanb;
      d = std::real(Md2(2,2)) + mYd(2)*mYd(2) + DdRR;
      /*
      cout<<"a = "<<a<<endl;
      cout<<"b = "<<b<<endl;
      cout<<"d = "<<d<<endl;
      */
      msbot_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
      msbot_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b) ;

      if (msbot_2 >=0.) msbot_2 = std::sqrt(msbot_2);
      else P = P*10.;

      if (msbot_1 >=0.) msbot_1 = std::sqrt(msbot_1);
      else P = P*10.;

      if (P!=1.) return P;

      s2sbot = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
      c2sbot = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);

      //ssbot2 = 0.5 - 0.5*c2sbot;
      //csbot2 = 1. - ssbot2;

      // check
      /*
      cout<<"mb_1 = "<<msbot_1<<"\n";
      cout<<"mb_2 = "<<msbot_2<<"\n";
      cout<<"std::sqrt(a*ctop2 + d*stop2 + b*s2top) = "<<std::sqrt(a*csbot2 + d*ssbot2 + b*s2sbot)<<"\n";
      cout<<"std::sqrt(a*stop2 + d*ctop2 - b*s2top) = "<<std::sqrt(a*ssbot2 + d*csbot2 - b*s2sbot)<<"\n";
      cout<<"b*c2top - 0.5*(a-d)*s2top = "<<b*c2sbot - 0.5*(a-d)*s2sbot<<"\n";
      cout<<"ctop2 = "<<csbot2<<"\n";
      cout<<"std::sqrt(ctop2) = "<<std::sqrt(csbot2)<<"\n";
      cout<<"stop2 = "<<ssbot2<<"\n";
      cout<<"std::sqrt(stop2) = "<<std::sqrt(ssbot2)<<"\n";
      */

      // double mstau_1, mstau_2, s2stau, c2stau, sstau2, cstau2;

      a = std::real(ML2(2,2)) + mYe(2)*mYe(2) + DeLL;
      b = vd*std::real(TYe(2,2)) - mu*mYe(2)*tanb;
      d = std::real(Me2(2,2)) + mYe(2)*mYe(2) + DeRR;
      /*
      cout<<"std::sqrt(a) = "<<std::sqrt(a)<<endl;
      cout<<"b = "<<b<<endl;
      cout<<"d = "<<std::sqrt(d)<<endl;
      */
      mstau_2 = 0.5*(a+d) + 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);
      mstau_1 = 0.5*(a+d) - 0.5*std::sqrt((a-d)*(a-d) + 4.*b*b);

      if (mstau_2 >=0.) mstau_2 = std::sqrt(mstau_2);
      else P = P*10.;

      if (mstau_1 >=0.) mstau_1 = std::sqrt(mstau_1);
      else P = P*10.;

      if (P!=1.) return P;

      s2stau = -2.*b/std::sqrt((a-d)*(a-d) + 4.*b*b);
      c2stau = -(a-d)/std::sqrt((a-d)*(a-d) + 4.*b*b);
      /*
      sstau2 = 0.5 - 0.5*c2stau;
      cstau2 = 1. - sstau2;

      // check

      cout<<"mb_1 = "<<mstau_1<<"\n";
      cout<<"mb_2 = "<<mstau_2<<"\n";
      cout<<"std::sqrt(a*ctop2 + d*stop2 + b*s2top) = "<<std::sqrt(a*cstau2 + d*sstau2 + b*s2stau)<<"\n";
      cout<<"std::sqrt(a*stop2 + d*ctop2 - b*s2top) = "<<std::sqrt(a*sstau2 + d*cstau2 - b*s2stau)<<"\n";
      cout<<"b*c2top - 0.5*(a-d)*s2top = "<<b*c2stau - 0.5*(a-d)*s2stau<<"\n";
      cout<<"ctop2 = "<<cstau2<<"\n";
      cout<<"std::sqrt(ctop2) = "<<std::sqrt(cstau2)<<"\n";
      cout<<"stop2 = "<<sstau2<<"\n";
      cout<<"std::sqrt(stop2) = "<<std::sqrt(sstau2)<<"\n";
      */

      // chargino mass matrix
      U(0, 0) = dcomplex(std::real(M(1)), 0.);
      U(0, 1) = g(1) * vu;
      U(1, 0) = g(1) * vd;
      U(1, 1) = mu;

      //cout<<"U =" << U;

      // MatrixXcd Ch_mass(2,2), Ch_diag(2,2);
      Ch_mass = U;

      bi_un_tr(U, 2, mMch, V);

      //  Matrix <double> U_mix(2,2), V_mix(2,2);
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
          U_mix(i,j) = std::real(U(i,j));
          V_mix(i,j) = std::real(V(j,i));
        }
      }
      V = (V).conjugate().eval();
      /*
      cout<<"mMch = "<<mMch<<endl;
      cout<<"U_mix = "<<U_mix<<endl;
      cout<<"V_mix = "<<V_mix<<endl;
      */

      // neutralino mass matrix
      Ntr(0, 0) = std::real(M[0]);
      Ntr(0, 1) = 0.;
      Ntr(0, 2) = - g[0]*std::sqrt(3./5) * vd/std::sqrt(2.);
      Ntr(0, 3) = g[0]*std::sqrt(3./5) * vu/std::sqrt(2.);

      Ntr(1, 1) = std::real(M[1]);
      Ntr(1, 2) = g[1] * vd/std::sqrt(2.);
      Ntr(1, 3) = - g[1] * vu/std::sqrt(2.);

      Ntr(2, 2) = 0.;
      Ntr(2, 3) = - mu;

      Ntr(3, 3) = 0.;

      for (int i = 1; i <= 3; i++) {
        for (int j = 0; j < i; j++) {
          Ntr(i, j) = Ntr(j, i);
        }
      }

      //cout<<"Ntr = "<<Ntr<<endl;

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          Ntr_zal(i,j) = Ntr(i,j);
          Ntr_mass(i,j) = Ntr(i,j);
        }
      }

      diag_r_s(Ntr, 4, mMnt);
      //cout<<"Ntr = " <<Ntr<<endl;
      //result is Ntr^T
      Nt_mix = Ntr;
      Ntr = (Ntr).transpose().eval();
      for (int i = 0; i < 4; i++) {
        if (mMnt(i) < 0.) {
          mMnt(i) = -mMnt(i);
          // forming P^*
          PNt(i) = dcomplex(0.,-1.);
        }
      }

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
          Nt(i,j) = PNt(i)*Ntr(i,j);
        }
      }

      Mg = std::real(M(2)) - ( - std::real(M(2))
           * (15. + 18.*log(MEW/std::real(M(2)))) )*g(2)*g(2)/(16.*pi*pi);

      for (int a = 0; a < 6; a++)
        Mg = Mg + std::real(M(2)) * ( B1(std::real(M(2)), 0., mMu(a)) +
             B1(std::real(M(2)), 0., mMd(a)) )*g(2)*g(2)/(16.*pi*pi);

      //  cout<<"Mg = "<<Mg<<endl;
    
    /*
   for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
        Yu(i,j) = Yu(i,j)*vu;
        Yd(i,j) = Yd(i,j)*vd;
        Ye(i,j) = Ye(i,j)*vd;
      }
    }
    */
    return P;
}
    //  MatrixXd G_stop(2,2), M_stop(2,2);

  
