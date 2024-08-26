#include "RGEsDefiner.h"

double KroneckerDelta(int i, int j) {
	return (i == j) ? 1.0 : 0.0;
}

// Defines the Beta functions -> dy(t)/dt and returns dydt as a vector:
VectorXcd RGEsystem(double t, const VectorXcd &y){
	double pi = M_PI;
	double LoopFactor = 1.0 / std::pow(4 * pi, 2.0);

	dcomplex MHu2, MHd2, Bmu;
	Matrix3cd Yu, Yd, Ye, TYu, TYd, TYe, Mq2, Ml2, Mu2, Md2, Me2;
	Vector3cd M;
    Vector3d g;
	//Storing the parameters into y(t) vector:
	// Mapping the parameters into a vector y(t):
	int k = 0;
	for(int i=0; i<=2; i++)
	{
		for (int j=0; j<=2; j++)
		{
			Yu(i,j) = y(k);
			Yd(i,j) = y(9+k);
			Ye(i,j) = y(18+k);
			TYu(i,j) = y(27+k);
			TYd(i,j) = y(36+k);
			TYe(i,j) = y(45+k);
			Mq2(i,j) = y(54+k);
			Ml2(i,j) = y(63+k);
			Mu2(i,j) = y(72+k);
			Md2(i,j) = y(81+k);
			Me2(i, j) = y(90+k);
			k++;
		}
		M(i) = y(99+i);
		g(i) = std::real(y(102+i));
	}
	MHu2 = y(105);
	MHd2 = y(106);
	//Bmu = y(107);
    
    Yu = Yu.transpose().eval();     // In model.cpp, Yukawas are defined with doublets on the left,
    Yd = Yd.transpose().eval();     // but RGEs are written with doublets on the right! We Transpose
    Ye = Ye.transpose().eval();     // them to be consistent!
                            
    TYu = TYu.transpose().eval();
    TYd = TYd.transpose().eval();
    TYe = TYe.transpose().eval();
    
    double Betag1l1, Betag1l2, Betag2l1, Betag2l2, Betag3l1, Betag3l2;
    
    dcomplex BetaM1l1, BetaM1l2, BetaM2l1, BetaM2l2, BetaM3l1, BetaM3l2, BetaMHd2l1, BetaMHd2l2, BetaMHu2l1, BetaMHu2l2;//BetaBmul1, BetaBmul2;
    
	Matrix3cd BetaYdl1, BetaYdl2, BetaYul1, BetaYul2, BetaYel1, BetaYel2, BetaTYdl1, BetaTYdl2, BetaTYul1, BetaTYul2,
					BetaTYel1, BetaTYel2, BetaMq2l1, BetaMq2l2, BetaMl2l1, BetaMl2l2, BetaMu2l1, BetaMu2l2, BetaMd2l1,
					BetaMd2l2, BetaMe2l1, BetaMe2l2;

	dcomplex dg1_dt, dg2_dt, dg3_dt, dM1_dt, dM2_dt, dM3_dt, dMHu2_dt, dMHd2_dt, dBmu_dt;
	Matrix3cd dYd_dt, dYu_dt, dYe_dt, dTYd_dt, dTYu_dt, dTYe_dt, dMq2_dt, dMl2_dt, dMu2_dt, dMd2_dt, dMe2_dt;


	Betag1l1 = (33.0*std::pow(g(0), 3.0))/5.0;
	Betag1l2 =  std::real((std::pow(g(0), 3.0)*(199.0*std::pow(g(0), 2.0)
					+ 135.0*std::pow(g(1), 2.0) + 440.0*std::pow(g(2), 2.0)
					- 70.0*(Yd *  Yd.adjoint()).trace() - 90.0*(Ye * Ye.adjoint()).trace()
					- 130.0*(Yu *  Yu.adjoint()).trace())) / 25.0);

	Betag2l1 = std::pow(g(1), 3.0);
	Betag2l2 = std::real((std::pow(g(1), 3.0)*(9.0*std::pow(g(0), 2.0) + 125.0*std::pow(g(1), 2.0)
					+ 120.0*std::pow(g(2), 2.0)
					- 30.0*(Yd * Yd.adjoint()).trace() - 10.0*(Ye * Ye.adjoint()).trace()
					- 30.0*(Yu * Yu.adjoint()).trace())) / 5.0);

	Betag3l1 = -3.0*std::pow(g(2), 3.0);
	Betag3l2 = std::real((std::pow(g(2), 3.0)*(11.0*std::pow(g(0), 2.0) + 45.0*std::pow(g(1), 2.0)
					+ 70.0*std::pow(g(2), 2.0)
					- 20.0*(Yd * Yd.adjoint()).trace() - 20.0*(Yu * Yu.adjoint()).trace())) / 5.0);
	BetaM1l1 = (66.0*std::pow(g(0), 2.0)*M(0))/5.0;
	BetaM1l2 = (2.0*std::pow(g(0), 2.0)*(398.0*std::pow(g(0), 2.0)*M(0)
					+ 135.0*std::pow(g(1), 2.0)*M(0)
					+ 440.0*std::pow(g(2), 2.0)*M(0)
					+ 440.0*std::pow(g(2), 2.0)*M(2) + 135.0*std::pow(g(1), 2.0)*M(1)
					- 70.0*M(0)* (Yd * Yd.adjoint() ).trace() - 90.0*M(0)* (Ye * Ye.adjoint()).trace()
					- 130.0*M(0)* (Yu * Yu.adjoint()).trace() + 70.0* (Yd.adjoint() * TYd).trace()
					+ 90.0* (Ye.adjoint() * TYe).trace() + 130.0* (Yu.adjoint() * TYu).trace()))/25.0;

	BetaM2l1 = 2.0*std::pow(g(1), 2.0)*M(1);
	BetaM2l2 = (2.0*std::pow(g(1), 2.0)*(9.0*std::pow(g(0), 2.0)*M(0)
					+ 120.0*std::pow(g(2), 2.0)*M(2)
					+ 9.0*std::pow(g(0), 2.0)*M(1) + 250.0*std::pow(g(1), 2.0)*M(1)
					+ 120.0*std::pow(g(2), 2.0)*M(1)
					- 30.0*M(1)* (Yd * Yd.adjoint()).trace() - 10.0*M(1)* (Ye * Ye.adjoint()).trace()
					- 30.0*M(1)* (Yu * Yu.adjoint()).trace() + 30.0* (Yd.adjoint() * TYd).trace()
					+ 10.0* (Ye.adjoint() * TYe).trace() + 30.0* (Yu.adjoint() * TYu).trace()))/5.0;

	BetaM3l1 = -6.0*std::pow(g(2), 2.0)*M(2);
	BetaM3l2 = (2.0*std::pow(g(2), 2.0)*(11.0*std::pow(g(0), 2.0)*M(0)
					+ 11.0*std::pow(g(0), 2.0)*M(2)
					+ 45.0*std::pow(g(1), 2.0)*M(2) + 140.0*std::pow(g(2), 2.0)*M(2)
					+ 45.0*std::pow(g(1), 2.0)*M(1) - 20.0*M(2)* (Yd * Yd.adjoint()).trace()
					- 20.0*M(2)* (Yu * Yu.adjoint()).trace() + 20.0* (Yd.adjoint() * TYd).trace()
					+ 20.0* (Yu.adjoint() * TYu).trace()))/5.0;

	BetaMHd2l1 = (-6.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0)))/5.0
					- 6.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1)) -
					std::sqrt(3.0/5.0)*g(0)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace()
					+ (Me2).trace() - (Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace()))
					+ 6.0*MHd2*(Yd * Yd.adjoint()).trace() +
					2.0*MHd2*(Ye * Ye.adjoint()).trace() + 6.0*(TYd.conjugate() * TYd.transpose()).trace() +
					2.0*(TYe.conjugate() * TYe.transpose()).trace() + 6.0*(Md2 * Yd * Yd.adjoint()).trace() +
					2.0*(Me2 * Ye * Ye.adjoint()).trace() + 2.0*(Ml2 * Ye.adjoint() * Ye).trace() +
					6.0*(Mq2 * Yd.adjoint() * Yd).trace();
	BetaMHd2l2 = (15.0*std::pow(g(1), 2.0)*(55.0*std::pow(g(1), 2.0)*M(1)
					+ 3.0*std::pow(g(0), 2.0)*(M(0) + 2.0*M(1)))*std::conj(M(1)) +
					std::pow(g(0), 2.0)*std::conj(M(0))*(621.0*std::pow(g(0), 2.0)*M(0)
					+ 90.0*std::pow(g(1), 2.0)*M(0) + 45.0*std::pow(g(1), 2.0)*M(1) -
					40.0*M(0)*(Yd * Yd.adjoint()).trace() + 120.0*M(0)*(Ye * Ye.adjoint()).trace() +
					20.0*(Yd.adjoint() * TYd).trace() - 60.0*(Ye.adjoint() * TYe).trace()) +
					10.0*(15.0*std::pow(g(1), 4.0)*((MHd2 + MHu2 + (Ml2).trace() + 3.0*(Mq2).trace())/2.0)
					+ 3.0*std::pow(g(0), 2.0)*((std::pow(g(0), 2.0)*(3.0*MHd2 + 3.0*MHu2
					+ 2.0*(Md2).trace() + 6.0*(Me2).trace() + 3.0*(Ml2).trace()
					+ (Mq2).trace() + 8.0*(Mu2).trace()))/10.0)
					- 2.0*std::sqrt(15.0)*g(0)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
					- 45.0*std::pow(g(1), 2.0)*MHd2 + 9.0*std::pow(g(0), 2.0)*MHu2
					+ 45.0*std::pow(g(1), 2.0)*MHu2 +
					4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Md2).trace()
					+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
					9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
					+ std::pow(g(0), 2.0)*(Mq2).trace() +
					45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
					- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
					160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
					30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
					60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
					- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
					60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
					+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
					120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
					- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
					(20.0*std::sqrt(15.0))) +
					(-2.0*std::pow(g(0), 2.0)*MHd2 + 80.0*std::pow(g(2), 2.0)*MHd2
					+ 160.0*std::pow(g(2), 2.0)*M(2)*std::conj(M(2)))*
					(Yd * Yd.adjoint()).trace() + 6.0*std::pow(g(0), 2.0)*MHd2*(Ye * Ye.adjoint()).trace() -
					80.0*std::pow(g(2), 2.0)*std::conj(M(2))*(Yd.adjoint() * TYd).trace() +
					2.0*std::pow(g(0), 2.0)*M(0)*(TYd.conjugate() * Yd.transpose()).trace()
					- 80.0*std::pow(g(2), 2.0)*M(2)*(TYd.conjugate() * Yd.transpose()).trace()
					- 2.0*std::pow(g(0), 2.0)*(TYd.conjugate() * TYd.transpose()).trace() +
					80.0*std::pow(g(2), 2.0)*(TYd.conjugate() * TYd.transpose()).trace()
					- 6.0*std::pow(g(0), 2.0)*M(0)*(TYe.conjugate() * Ye.transpose()).trace()
					+ 6.0*std::pow(g(0), 2.0)*(TYe.conjugate() * TYe.transpose()).trace() -
					2.0*std::pow(g(0), 2.0)*(Md2 * Yd * Yd.adjoint()).trace()
					+ 80.0*std::pow(g(2), 2.0)*(Md2 * Yd * Yd.adjoint()).trace() +
					6.0*std::pow(g(0), 2.0)*(Me2 * Ye * Ye.adjoint()).trace()
					+ 6.0*std::pow(g(0), 2.0)*(Ml2 * Ye.adjoint() * Ye).trace() -
					2.0*std::pow(g(0), 2.0)*(Mq2 * Yd.adjoint() * Yd).trace()
					+ 80.0*std::pow(g(2), 2.0)*(Mq2 * Yd.adjoint() * Yd).trace() -
					90.0*MHd2*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace() -
					90.0*(Yd * Yd.adjoint() * TYd * TYd.adjoint()).trace() -
					15.0*MHd2*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					15.0*MHu2*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					15.0*(Yd * Yu.adjoint() * TYu * TYd.adjoint()).trace() -
					90.0*(Yd * TYd.adjoint() * TYd * Yd.adjoint()).trace() -
					15.0*(Yd * TYu.adjoint() * TYu * Yd.adjoint()).trace() -
					30.0*MHd2*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace() -
					30.0*(Ye * Ye.adjoint() * TYe * TYe.adjoint()).trace() -
					30.0*(Ye * TYe.adjoint() * TYe * Ye.adjoint()).trace() -
					15.0*(Yu * Yd.adjoint() * TYd * TYu.adjoint()).trace() -
					15.0*(Yu * TYd.adjoint() * TYd * Yu.adjoint()).trace() -
					90.0*(Md2 * Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace() -
					15.0*(Md2 * Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					30.0*(Me2 * Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace() -
					30.0*(Ml2 * Ye.adjoint() * Ye * Ye.adjoint() * Ye).trace() -
					90.0*(Mq2 * Yd.adjoint() * Yd * Yd.adjoint() * Yd).trace() -
					15.0*(Mq2 * Yd.adjoint() * Yd * Yu.adjoint() * Yu).trace() -
					15.0*(Mq2 * Yu.adjoint() * Yu * Yd.adjoint() * Yd).trace() -
					15.0*(Mu2 * Yu * Yd.adjoint() * Yd * Yu.adjoint()).trace()))/25.0;

	BetaMHu2l1 = (-6.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0)))/5.0
					- 6.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1)) +
					std::sqrt(3.0/5.0)*g(0)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace() + (Me2).trace() -
					(Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace())) + 6.0*MHu2*(Yu * Yu.adjoint()).trace() +
					6.0*(TYu.conjugate() * TYu.transpose()).trace() + 6.0*(Mq2 * Yu.adjoint() * Yu).trace() +
					6.0*(Mu2 * Yu * Yu.adjoint()).trace();
	BetaMHu2l2 = 33.0*std::pow(g(1), 4.0)*M(1)*std::conj(M(1)) +
					(9.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*(M(0)
					+ 2.0*M(1))*std::conj(M(1)))/5.0
					+ 6.0*std::pow(g(1), 4.0)*((MHd2 + MHu2 + (Ml2).trace()
					+ 3.0*(Mq2).trace())/2.0) + (6.0*std::pow(g(0), 2.0)*((std::pow(g(0), 2.0)*(3.0*MHd2
					+ 3.0*MHu2 + 2.0*(Md2).trace() + 6.0*(Me2).trace() +
					3.0*(Ml2).trace() + (Mq2).trace() + 8.0*(Mu2).trace()))/10.0))/5.0
					+ 4.0*std::sqrt(3.0/5.0)*g(0)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
					- 45.0*std::pow(g(1), 2.0)*MHd2
					+ 9.0*std::pow(g(0), 2.0)*MHu2 + 45.0*std::pow(g(1), 2.0)*MHu2 +
					4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Md2).trace()
					+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
					9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
					+ std::pow(g(0), 2.0)*(Mq2).trace() +
					45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
					- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
					160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
					30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
					60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
					- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
					60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
					+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
					120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
					- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
					(20.0*std::sqrt(15.0))) +
					(8.0*std::pow(g(0), 2.0)*MHu2*(Yu * Yu.adjoint()).trace())/5.0
					+ 32.0*std::pow(g(2), 2.0)*MHu2*(Yu * Yu.adjoint()).trace() +
					64.0*std::pow(g(2), 2.0)*M(2)*std::conj(M(2))*(Yu * Yu.adjoint()).trace() +
					(std::pow(g(0), 2.0)*std::conj(M(0))*(621.0*std::pow(g(0), 2.0)*M(0)
					+ 90.0*std::pow(g(1), 2.0)*M(0) + 45.0*std::pow(g(1), 2.0)*M(1) +
					80.0*M(0)*(Yu * Yu.adjoint()).trace() - 40.0*(Yu.adjoint() * TYu).trace()))/25.0 -
					32.0*std::pow(g(2), 2.0)*std::conj(M(2))*(Yu.adjoint() * TYu).trace() -
					(8.0*std::pow(g(0), 2.0)*M(0)*(TYu.conjugate() * Yu.transpose()).trace())/5.0 -
					32.0*std::pow(g(2), 2.0)*M(2)*(TYu.conjugate() * Yu.transpose()).trace() +
					(8.0*std::pow(g(0), 2.0)*(TYu.conjugate() * TYu.transpose()).trace())/5.0 +
					32.0*std::pow(g(2), 2.0)*(TYu.conjugate() * TYu.transpose()).trace()
					+ (8.0*std::pow(g(0), 2.0)*(Mq2 * Yu.adjoint() * Yu).trace())/
					5.0 + 32.0*std::pow(g(2), 2.0)*(Mq2 * Yu.adjoint() * Yu).trace()
					+ (8.0*std::pow(g(0), 2.0)*(Mu2 * Yu * Yu.adjoint()).trace())/
					5.0 + 32.0*std::pow(g(2), 2.0)*(Mu2 * Yu * Yu.adjoint()).trace() -
					6.0*MHd2*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					6.0*MHu2*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					6.0*(Yd * Yu.adjoint() * TYu * TYd.adjoint()).trace() -
					6.0*(Yd * TYu.adjoint() * TYu * Yd.adjoint()).trace() -
					6.0*(Yu * Yd.adjoint() * TYd * TYu.adjoint()).trace() -
					36.0*MHu2*(Yu * Yu.adjoint() * Yu * Yu.adjoint()).trace() -
					36.0*(Yu * Yu.adjoint() * TYu * TYu.adjoint()).trace() -
					6.0*(Yu * TYd.adjoint() * TYd * Yu.adjoint()).trace() -
					36.0*(Yu * TYu.adjoint() * TYu * Yu.adjoint()).trace() -
					6.0*(Md2 * Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					6.0*(Mq2 * Yd.adjoint() * Yd * Yu.adjoint() * Yu).trace() -
					6.0*(Mq2 * Yu.adjoint() * Yu * Yd.adjoint() * Yd).trace() -
					36.0*(Mq2 * Yu.adjoint() * Yu * Yu.adjoint() * Yu).trace() -
					6.0*(Mu2 * Yu * Yd.adjoint() * Yd * Yu.adjoint()).trace() -
					36.0*(Mu2 * Yu * Yu.adjoint() * Yu * Yu.adjoint()).trace();

	/*BetaBmul1 = (6.0*std::pow(g(0), 2.0)*M(0)*mu)/5.0 + 6.0*std::pow(g(1), 2.0)*M(1)*mu +
					Bmu*((-3.0*std::pow(g(0), 2.0))/5.0 - 3.0*std::pow(g(1), 2.0)
					+ 3.0*(Yd * Yd.adjoint()).trace() +
					(Ye * Ye.adjoint()).trace() + 3.0*(Yu * Yu.adjoint()).trace()) +
					6.0*mu*(Yd.adjoint() * TYd).trace() + 2.0*mu*(Ye.adjoint() * TYe).trace() +
					6.0*mu*(Yu.adjoint() * TYu).trace();
	BetaBmul2 = Bmu *((207.0*std::pow(g(0), 4.0))/50.0
					+ (9.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0))/5.0
					+ (15.0*std::pow(g(1), 4.0))/2.0 -
					(2.0*(std::pow(g(0), 2.0) - 40.0*std::pow(g(2), 2.0))*(Yd * Yd.adjoint()).trace())/5.0
					+ (6.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint()).trace())/
					5.0 + (4.0*std::pow(g(0), 2.0)*(Yu * Yu.adjoint()).trace())/5.0
					+ 16.0*std::pow(g(2), 2.0)*(Yu * Yu.adjoint()).trace() -
					9.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace()
					- 6.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace() -
					3.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace()
					- 9.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint()).trace()) -
					(2.0*mu*(207.0*std::pow(g(0), 4.0)*M(0)
					+ 45.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(0)
					+ 45.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(1) +
					375.0*std::pow(g(1), 4.0)*M(1) - 10.0*(std::pow(g(0), 2.0)*M(0)
					- 40.0*std::pow(g(2), 2.0)*M(2))*(Yd * Yd.adjoint()).trace() +
					30.0*std::pow(g(0), 2.0)*M(0)*(Ye * Ye.adjoint()).trace()
					+ 20.0*std::pow(g(0), 2.0)*M(0)*(Yu * Yu.adjoint()).trace() +
					400.0*std::pow(g(2), 2.0)*M(2)*(Yu * Yu.adjoint()).trace()
					+ 10.0*std::pow(g(0), 2.0)*(Yd.adjoint() * TYd).trace() -
					400.0*std::pow(g(2), 2.0)*(Yd.adjoint() * TYd).trace()
					- 30.0*std::pow(g(0), 2.0)*(Ye.adjoint() * TYe).trace() -
					20.0*std::pow(g(0), 2.0)*(Yu.adjoint() * TYu).trace()
					- 400.0*std::pow(g(2), 2.0)*(Yu.adjoint() * TYu).trace() +
					450.0*(Yd * Yd.adjoint() * TYd * Yd.adjoint()).trace() +
					150.0*(Yd * Yu.adjoint() * TYu * Yd.adjoint()).trace() +
					150.0*(Ye * Ye.adjoint() * TYe * Ye.adjoint()).trace() +
					150.0*(Yu * Yd.adjoint() * TYd * Yu.adjoint()).trace() +
					450.0*(Yu * Yu.adjoint() * TYu * Yu.adjoint()).trace()))/25.0;*/
	k = 0;
	for(int i=0; i<=2; i++)
	{
		for (int j=0.0; j<=2.0; j++) {
			BetaYdl1(i, j) = ((-7.0*std::pow(g(0), 2.0))/15.0 - 3.0*std::pow(g(1), 2.0)
							- (16.0*std::pow(g(2), 2.0))/3.0 + 3.0*(Yd * Yd.adjoint()).trace()
							+ (Ye * Ye.adjoint()).trace())*Yd(i, j) + 3.0*(Yd * Yd.adjoint() * Yd)(i, j)
							+ (Yd * Yu.adjoint() * Yu)(i, j);
			BetaYdl2(i, j) = ((287.0*std::pow(g(0), 4.0))/90.0
							+ std::pow(g(0), 2.0)*std::pow(g(1), 2.0)
							+ (15.0*std::pow(g(1), 4.0))/2.0
							+ (8.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0))/9.0
							+ 8.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)
							- (16.0*std::pow(g(2), 4.0))/9.0 - (2.0*(std::pow(g(0), 2.0)
							- 40.0*std::pow(g(2), 2.0))*(Yd * Yd.adjoint()).trace())/5.0
							+ (6.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint()).trace())/5.0
							- 9.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace()
							- 3.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()
							- 3.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace()) * Yd(i, j)
							+ ((4.0*std::pow(g(0), 2.0))/5.0 + 6.0*std::pow(g(1), 2.0)
							- 9.0*(Yd * Yd.adjoint()).trace()
							- 3.0*(Ye * Ye.adjoint()).trace())*(Yd * Yd.adjoint() * Yd)(i, j)
							+ (4.0*std::pow(g(0), 2.0)*(Yd * Yu.adjoint() * Yu)(i, j))/5.0
							- 3.0*(Yu * Yu.adjoint()).trace()*(Yd * Yu.adjoint() * Yu)(i, j)
							- 4.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j)
							- 2.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint() * Yd)(i, j)
							- 2.0*(Yd * Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j);

			BetaYel1(i, j) = ((-9.0*std::pow(g(0), 2.0))/5.0 - 3.0*std::pow(g(1), 2.0)
							+ 3.0*(Yd * Yd.adjoint()).trace()
							+ (Ye * Ye.adjoint()).trace())*Ye(i, j) + 3.0*(Ye * Ye.adjoint() * Ye)(i, j);
			BetaYel2(i, j) = ((-4.0*(std::pow(g(0), 2.0) - 40.0*std::pow(g(2), 2.0))*(Yd * Yd.adjoint()).trace()
							+ 3.0*(45.0*std::pow(g(0), 4.0) + 6.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)
							+ 25.0*std::pow(g(1), 4.0) + 4.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint()).trace()
							- 30.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace()
							- 10.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()
							- 10.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace()))*Ye(i, j))/10.0
							+ (6.0*std::pow(g(1), 2.0) - 9.0*(Yd * Yd.adjoint()).trace()
							- 3.0*(Ye * Ye.adjoint()).trace()) * (Ye * Ye.adjoint() * Ye)(i, j)
							- 4.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint() * Ye)(i, j);

			BetaYul1(i, j) = -1.0/15.0*((13.0*std::pow(g(0), 2.0)
							+ 45.0*std::pow(g(1), 2.0) + 80.0*std::pow(g(2), 2.0)
							- 45.0*(Yu * Yu.adjoint()).trace())*Yu(i, j))
							+ (Yu * Yd.adjoint() * Yd)(i, j) + 3.0*(Yu * Yu.adjoint() * Yu)(i, j);
			BetaYul2(i, j) = ((2743.0*std::pow(g(0), 4.0))/450.0
			  				+ std::pow(g(0), 2.0)*std::pow(g(1), 2.0)
							+ (15.0*std::pow(g(1), 4.0))/2.0 + (136.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0))/45.0
							+ 8.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0) - (16.0*std::pow(g(2), 4.0))/9.0
							+ (4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Yu * Yu.adjoint()).trace())/5.0
							- 3.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()
							- 9.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint()).trace()) * Yu(i, j)
							+ ((2.0*std::pow(g(0), 2.0))/5.0 - 3.0*(Yd * Yd.adjoint()).trace()
							- (Ye * Ye.adjoint()).trace()) * (Yu * Yd.adjoint() * Yd)(i, j)
							+ (2.0*std::pow(g(0), 2.0)*(Yu * Yu.adjoint() * Yu)(i, j))/5.0
							+ 6.0*std::pow(g(1), 2.0)*(Yu * Yu.adjoint() * Yu)(i, j)
							- 9.0*(Yu * Yu.adjoint()).trace() * (Yu * Yu.adjoint() * Yu)(i, j)
							- 2.0*(Yu * Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j)
							- 2.0*(Yu * Yd.adjoint() * Yd * Yu.adjoint() * Yu)(i, j)
							- 4.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j);

			BetaTYdl1(i, j) = ((14.0*std::pow(g(0), 2.0)*M(0))/15.0 + (32.0*std::pow(g(2), 2.0)*M(2))/3.0
							+ 6.0*std::pow(g(1), 2.0)*M(1) + 6.0*(Yd.adjoint() * TYd).trace()
							+ 2.0*(Ye.adjoint() * TYe).trace())*Yd(i, j)
							+ 4.0*(Yd * Yd.adjoint() * TYd)(i, j) + 2.0*(Yd * Yu.adjoint() * TYu)(i, j)
							+ 5.0*(TYd * Yd.adjoint() * Yd)(i, j) + (TYd * Yu.adjoint() * Yu)(i, j)
							- (7.0*std::pow(g(0), 2.0)*TYd(i, j))/15.0 - 3.0*std::pow(g(1), 2.0)*TYd(i, j)
							- (16.0*std::pow(g(2), 2.0)*TYd(i, j))/3.0
							+ 3.0*(Yd * Yd.adjoint()).trace()*TYd(i, j)
							+ (Ye * Ye.adjoint()).trace()*TYd(i, j);
			BetaTYdl2(i, j) = (-2.0*(287.0*std::pow(g(0), 4.0)*M(0) + 45.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(0)
							+ 40.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*M(0)
							+ 40.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*M(2)
							+ 360.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(2)
							- 160.0*std::pow(g(2), 4.0)*M(2) + 45.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(1)
							+ 675.0*std::pow(g(1), 4.0)*M(1) + 360.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(1)
							- 18.0*(std::pow(g(0), 2.0)*M(0)
							- 40.0*std::pow(g(2), 2.0)*M(2))*(Yd * Yd.adjoint()).trace()
							+ 54.0*std::pow(g(0), 2.0)*M(0)*(Ye * Ye.adjoint()).trace()
							+ 18.0*std::pow(g(0), 2.0)*(Yd.adjoint() * TYd).trace()
							- 720.0*std::pow(g(2), 2.0)*(Yd.adjoint() * TYd).trace()
							- 54.0*std::pow(g(0), 2.0)*(Ye.adjoint() * TYe).trace()
							+ 810.0*(Yd * Yd.adjoint() * TYd * Yd.adjoint()).trace()
							+ 135.0*(Yd * Yu.adjoint() * TYu * Yd.adjoint()).trace()
							+ 270.0*(Ye * Ye.adjoint() * TYe * Ye.adjoint()).trace()
							+ 135.0*(Yu * Yd.adjoint() * TYd * Yu.adjoint()).trace())*Yd(i, j))/45.0
							- (2.0*(4.0*std::pow(g(0), 2.0)*M(0)
							+ 30.0*std::pow(g(1), 2.0)*M(1) + 45.0*(Yd.adjoint() * TYd).trace()
							+ 15.0*(Ye.adjoint() * TYe).trace())*(Yd * Yd.adjoint() * Yd)(i, j))/5.0
							+ (6.0*std::pow(g(0), 2.0)*(Yd * Yd.adjoint() * TYd)(i, j))/5.0
							+ 6.0*std::pow(g(1), 2.0)*(Yd * Yd.adjoint() * TYd)(i, j)
							- 12.0*(Yd * Yd.adjoint()).trace()*(Yd * Yd.adjoint() * TYd)(i, j)
							- 4.0*(Ye * Ye.adjoint()).trace()*(Yd * Yd.adjoint() * TYd)(i, j)
							- (8.0*std::pow(g(0), 2.0)*M(0)*(Yd * Yu.adjoint() * Yu)(i, j))/5.0
							- 6.0*(Yu.adjoint() * TYu).trace()*(Yd * Yu.adjoint() * Yu)(i, j)
							+ (8.0*std::pow(g(0), 2.0)*(Yd * Yu.adjoint() * TYu)(i, j))/5.0
							- 6.0*(Yu * Yu.adjoint()).trace()*(Yd * Yu.adjoint() * TYu)(i, j)
							+ (6.0*std::pow(g(0), 2.0)*(TYd * Yd.adjoint() * Yd)(i, j))/5.0
							+ 12.0*std::pow(g(1), 2.0)*(TYd * Yd.adjoint() * Yd)(i, j)
							- 15.0*(Yd * Yd.adjoint()).trace()*(TYd * Yd.adjoint() * Yd)(i, j)
							- 5.0*(Ye * Ye.adjoint()).trace()*(TYd * Yd.adjoint() * Yd)(i, j)
							+ (4.0*std::pow(g(0), 2.0)*(TYd * Yu.adjoint() * Yu)(i, j))/5.0
							- 3.0*(Yu * Yu.adjoint()).trace()*(TYd * Yu.adjoint() * Yu)(i, j)
							- 6.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint() * TYd)(i, j)
							- 8.0*(Yd * Yd.adjoint() * TYd * Yd.adjoint() * Yd)(i, j)
							- 2.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint() * TYd)(i, j)
							- 4.0*(Yd * Yu.adjoint() * Yu * Yu.adjoint() * TYu)(i, j)
							- 4.0*(Yd * Yu.adjoint() * TYu * Yd.adjoint() * Yd)(i, j)
							- 4.0*(Yd * Yu.adjoint() * TYu * Yu.adjoint() * Yu)(i, j)
							- 6.0*(TYd * Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j)
							- 4.0*(TYd * Yu.adjoint() * Yu * Yd.adjoint() * Yd)(i, j)
							- 2.0*(TYd * Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j)
							+ (287.0*std::pow(g(0), 4.0)*TYd(i, j))/90.0
							+ std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*TYd(i, j)
							+ (15.0*std::pow(g(1), 4.0)*TYd(i, j))/2.0
							+ (8.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*TYd(i, j))/9.0
							+ 8.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*TYd(i, j)
							- (16.0*std::pow(g(2), 4.0)*TYd(i, j))/9.0
							- (2.0*std::pow(g(0), 2.0)*(Yd * Yd.adjoint()).trace()*TYd(i, j))/5.0
							+ 16.0*std::pow(g(2), 2.0)*(Yd * Yd.adjoint()).trace()*TYd(i, j)
							+ (6.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint()).trace()*TYd(i, j))/5.0
							- 9.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace()*TYd(i, j)
							- 3.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()*TYd(i, j)
							- 3.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace()*TYd(i, j);

			BetaTYel1(i, j) = ((18.0*std::pow(g(0), 2.0)*M(0))/5.0 + 6.0*std::pow(g(1), 2.0)*M(1)
							+ 6.0*(Yd.adjoint() * TYd).trace() + 2.0*(Ye.adjoint() * TYe).trace())*Ye(i, j)
							+ 4.0*(Ye * Ye.adjoint() * TYe)(i, j)
							+ 5.0*(TYe * Ye.adjoint() * Ye)(i, j)
							- (9.0*std::pow(g(0), 2.0)*TYe(i, j))/5.0
							- 3.0*std::pow(g(1), 2.0)*TYe(i, j) + 3.0*(Yd * Yd.adjoint()).trace()*TYe(i, j)
							+ (Ye * Ye.adjoint()).trace()*TYe(i, j);
			BetaTYel2(i, j) = (-2.0*(135.0*std::pow(g(0), 4.0)*M(0)
                            + 9.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(0)
							+ 9.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(1)
							+ 75.0*std::pow(g(1), 4.0)*M(1) + (-2.0*std::pow(g(0), 2.0)*M(0)
							+ 80.0*std::pow(g(2), 2.0)*M(2))*(Yd * Yd.adjoint()).trace()
							+ 6.0*std::pow(g(0), 2.0)*M(0)*(Ye * Ye.adjoint()).trace()
							+ 2.0*std::pow(g(0), 2.0)*(Yd.adjoint() * TYd).trace()
							- 80.0*std::pow(g(2), 2.0)*(Yd.adjoint() * TYd).trace()
							- 6.0*std::pow(g(0), 2.0)*(Ye.adjoint() * TYe).trace()
							+ 90.0*(Yd * Yd.adjoint() * TYd * Yd.adjoint()).trace()
							+ 15.0*(Yd * Yu.adjoint() * TYu * Yd.adjoint()).trace()
							+ 30.0*(Ye * Ye.adjoint() * TYe * Ye.adjoint()).trace()
							+ 15.0*(Yu * Yd.adjoint() * TYd * Yu.adjoint()).trace())*Ye(i, j))/5.0
							- 6.0*(2.0*std::pow(g(1), 2.0)*M(1) + 3.0*(Yd.adjoint() * TYd).trace()
							+ (Ye.adjoint() * TYe).trace())*(Ye * Ye.adjoint() * Ye)(i, j)
							+ (6.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint() * TYe)(i, j))/5.0
							+ 6.0*std::pow(g(1), 2.0)*(Ye * Ye.adjoint() * TYe)(i, j)
							- 12.0*(Yd * Yd.adjoint()).trace()*(Ye * Ye.adjoint() * TYe)(i, j)
							- 4.0*(Ye * Ye.adjoint()).trace()*(Ye * Ye.adjoint() * TYe)(i, j)
							- (6.0*std::pow(g(0), 2.0)*(TYe * Ye.adjoint() * Ye)(i, j))/5.0
							+ 12.0*std::pow(g(1), 2.0)*(TYe * Ye.adjoint() * Ye)(i, j)
							- 15.0*(Yd * Yd.adjoint()).trace()*(TYe * Ye.adjoint() * Ye)(i, j)
							- 5.0*(Ye * Ye.adjoint()).trace()*(TYe * Ye.adjoint() * Ye)(i, j)
							- 6.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint() * TYe)(i, j)
							- 8.0*(Ye * Ye.adjoint() * TYe * Ye.adjoint() * Ye)(i, j)
							- 6.0*(TYe * Ye.adjoint() * Ye * Ye.adjoint() * Ye)(i, j)
							+ (27.0*std::pow(g(0), 4.0)*TYe(i, j))/2.0
							+ (9.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*TYe(i, j))/5.0
							+ (15.0*std::pow(g(1), 4.0)*TYe(i, j))/2.0
							- (2.0*std::pow(g(0), 2.0)*(Yd * Yd.adjoint()).trace()*TYe(i, j))/5.0
							+ 16.0*std::pow(g(2), 2.0)*(Yd * Yd.adjoint()).trace()*TYe(i, j)
							+ (6.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint()).trace()*TYe(i, j))/5.0
							- 9.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint()).trace()*TYe(i, j)
							- 3.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()*TYe(i, j)
							- 3.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint()).trace()*TYe(i, j);

			BetaTYul1(i, j) = ((26.0*std::pow(g(0), 2.0)*M(0))/15.0 + (32.0*std::pow(g(2), 2.0)*M(2))/3.0
							+ 6.0*std::pow(g(1), 2.0)*M(1) + 6.0*(Yu.adjoint() * TYu).trace())*Yu(i, j)
							+ 2.0*(Yu * Yd.adjoint() * TYd)(i, j) + 4.0*(Yu * Yu.adjoint() * TYu)(i, j)
							+ (TYu * Yd.adjoint() * Yd)(i, j) + 5.0*(TYu * Yu.adjoint() * Yu)(i, j)
							- (13.0*std::pow(g(0), 2.0)*TYu(i, j))/15.0 - 3.0*std::pow(g(1), 2.0)*TYu(i, j)
							- (16.0*std::pow(g(2), 2.0)*TYu(i, j))/3.0 + 3.0*(Yu * Yu.adjoint()).trace()*TYu(i, j);
			BetaTYul2(i, j) = (-2.0*(2743.0*std::pow(g(0), 4.0)*M(0) + 225.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(0)
							+ 680.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*M(0)
							+ 680.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*M(2)
							+ 1800.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(2) - 800.0*std::pow(g(2), 4.0)*M(2)
							+ 225.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(1)
							+ 3375.0*std::pow(g(1), 4.0)*M(1) + 1800.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(1)
							+ 180.0*(std::pow(g(0), 2.0)*M(0) + 20.0*std::pow(g(2), 2.0)*M(2))*(Yu * Yu.adjoint()).trace()
							- 180.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Yu.adjoint() * TYu).trace()
							+ 675.0*(Yd * Yu.adjoint() * TYu * Yd.adjoint()).trace()
							+ 675.0*(Yu * Yd.adjoint() * TYd * Yu.adjoint()).trace()
							+ 4050.0*(Yu * Yu.adjoint() * TYu * Yu.adjoint()).trace())*Yu(i, j))/225.0
							- (2.0*(2.0*std::pow(g(0), 2.0)*M(0) + 15.0*(Yd.adjoint() * TYd).trace()
							+ 5.0*(Ye.adjoint() * TYe).trace())*(Yu * Yd.adjoint() * Yd)(i, j))/5.0
							+ (4.0*std::pow(g(0), 2.0)*(Yu * Yd.adjoint() * TYd)(i, j))/5.0
							- 6.0*(Yd * Yd.adjoint()).trace()*(Yu * Yd.adjoint() * TYd)(i, j)
							- 2.0*(Ye * Ye.adjoint()).trace()*(Yu * Yd.adjoint() * TYd)(i, j)
							- (4.0*std::pow(g(0), 2.0)*M(0)*(Yu * Yu.adjoint() * Yu)(i, j))/5.0
							- 12.0*std::pow(g(1), 2.0)*M(1)*(Yu * Yu.adjoint() * Yu)(i, j)
							- 18.0*(Yu.adjoint() * TYu).trace()*(Yu * Yu.adjoint() * Yu)(i, j)
							+ (6.0*std::pow(g(0), 2.0)*(Yu * Yu.adjoint() * TYu)(i, j))/5.0
							+ 6.0*std::pow(g(1), 2.0)*(Yu * Yu.adjoint() * TYu)(i, j)
							- 12.0*(Yu * Yu.adjoint()).trace()*(Yu * Yu.adjoint() * TYu)(i, j)
							+ (2.0*std::pow(g(0), 2.0)*(TYu * Yd.adjoint() * Yd)(i, j))/5.0
							- 3.0*(Yd * Yd.adjoint()).trace()*(TYu * Yd.adjoint() * Yd)(i, j)
							- (Ye * Ye.adjoint()).trace()*(TYu * Yd.adjoint() * Yd)(i, j)
							+ 12.0*std::pow(g(1), 2.0)*(TYu * Yu.adjoint() * Yu)(i, j)
							- 15.0*(Yu * Yu.adjoint()).trace()*(TYu * Yu.adjoint() * Yu)(i, j)
							- 4.0*(Yu * Yd.adjoint() * Yd * Yd.adjoint() * TYd)(i, j)
							- 2.0*(Yu * Yd.adjoint() * Yd * Yu.adjoint() * TYu)(i, j)
							- 4.0*(Yu * Yd.adjoint() * TYd * Yd.adjoint() * Yd)(i, j)
							- 4.0*(Yu * Yd.adjoint() * TYd * Yu.adjoint() * Yu)(i, j)
							- 6.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint() * TYu)(i, j)
							- 8.0*(Yu * Yu.adjoint() * TYu * Yu.adjoint() * Yu)(i, j)
							- 2.0*(TYu * Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j)
							- 4.0*(TYu * Yd.adjoint() * Yd * Yu.adjoint() * Yu)(i, j)
							- 6.0*(TYu * Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j)
							+ (2743.0*std::pow(g(0), 4.0)*TYu(i, j))/450.0
							+ std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*TYu(i, j)
							+ (15.0*std::pow(g(1), 4.0)*TYu(i, j))/2.0
							+ (136.0*std::pow(g(0), 2.0)*std::pow(g(2), 2.0)*TYu(i, j))/45.0
							+ 8.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*TYu(i, j)
							- (16.0*std::pow(g(2), 4.0)*TYu(i, j))/9.0
							+ (4.0*std::pow(g(0), 2.0)*(Yu * Yu.adjoint()).trace()*TYu(i, j))/5.0
							+ 16.0*std::pow(g(2), 2.0)*(Yu * Yu.adjoint()).trace()*TYu(i, j)
							- 3.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint()).trace()*TYu(i, j)
							- 9.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint()).trace()*TYu(i, j);

			BetaMq2l1(i, j) = (-2.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0))*KroneckerDelta(i, j))/15.0 -
							(32.0*std::pow(g(2), 2.0)*M(2)*std::conj(M(2))*KroneckerDelta(i, j))/3.0 -
							6.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))*KroneckerDelta(i, j) +
							(g(0)*KroneckerDelta(i, j)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace()
							+ (Me2).trace() - (Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace())))/std::sqrt(15.0) +
							2.0*MHd2*(Yd.adjoint() * Yd)(i, j) + 2.0*MHu2*(Yu.adjoint() * Yu)(i, j) +
							2.0*(TYd.adjoint() * TYd)(i, j) +
							2.0*(TYu.adjoint() * TYu)(i, j) + (Mq2 * Yd.adjoint() * Yd)(i, j) +
							(Mq2 * Yu.adjoint() * Yu)(i, j) + 2.0*(Yd.adjoint() * Md2 * Yd)(i, j) +
							(Yd.adjoint() * Yd * Mq2)(i, j) + 2.0*(Yu.adjoint() * Mu2 * Yu)(i, j) +
							(Yu.adjoint() * Yu * Mq2)(i, j);
			BetaMq2l2(i, j) = (-16.0*std::pow(g(2), 2.0)*(120.0*std::pow(g(2), 2.0)*M(2)
							- std::pow(g(0), 2.0)*(M(0) + 2.0*M(2)) -
							45.0*std::pow(g(1), 2.0)*(2.0*M(2)
							+ M(1)))*std::conj(M(2))*KroneckerDelta(i, j))/45.0 +
							(std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(0)*std::conj(M(1))*KroneckerDelta(i, j))/5.0 +
							16.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(2)*std::conj(M(1))*KroneckerDelta(i, j) +
							(2.0*std::pow(g(0), 2.0)*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))*KroneckerDelta(i, j))/5.0 +
							33.0*std::pow(g(1), 4.0)*M(1)*std::conj(M(1))*KroneckerDelta(i, j) +
							32.0*std::pow(g(1), 2.0)*std::pow(g(2), 2.0)*M(1)*std::conj(M(1))*KroneckerDelta(i, j) +
							6.0*std::pow(g(1), 4.0)*KroneckerDelta(i, j)*((MHd2 + MHu2 + (Ml2).trace()
							+ 3.0*(Mq2).trace())/2.0) + (32.0*std::pow(g(2), 4.0)*KroneckerDelta(i, j)*(((Md2).trace()
							+ 2.0*(Mq2).trace() + (Mu2).trace())/2.0))/3.0 +
							(2.0*std::pow(g(0), 2.0)*KroneckerDelta(i, j)*((std::pow(g(0), 2.0)*(3.0*MHd2
							+ 3.0*MHu2 + 2.0*(Md2).trace() + 6.0*(Me2).trace() +
							3.0*(Ml2).trace() + (Mq2).trace() + 8.0*(Mu2).trace()))/10.0))/15.0 +
							(4.0*g(0) * KroneckerDelta(i, j)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
							- 45.0*std::pow(g(1), 2.0)*MHd2 + 9.0*std::pow(g(0), 2.0)*MHu2
							+ 45.0*std::pow(g(1), 2.0)*MHu2 +
							4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Md2).trace()
							+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
							9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
							+ std::pow(g(0), 2.0)*(Mq2).trace() +
							45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
							- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
							160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
							30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
							60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
							- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
							60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
							+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
							120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
							- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
							(20.0*std::sqrt(15.0))))/std::sqrt(15.0) +
							(4.0*std::pow(g(0), 2.0)*MHd2*(Yd.adjoint() * Yd)(i, j))/5.0
							- 12.0*MHd2*(Yd * Yd.adjoint()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 4.0*MHd2*(Ye * Ye.adjoint()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 6.0*(TYd.conjugate() * TYd.transpose()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 2.0*(TYe.conjugate() * TYe.transpose()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 6.0*(Md2 * Yd * Yd.adjoint()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 2.0*(Me2 * Ye * Ye.adjoint()).trace()*
							(Yd.adjoint() * Yd)(i, j) - 2.0*(Ml2 * Ye.adjoint() * Ye).trace()*
							(Yd.adjoint() * Yd)(i, j) - 6.0*(Mq2 * Yd.adjoint() * Yd).trace()*
							(Yd.adjoint() * Yd)(i, j) - 6.0*(TYd.conjugate() * Yd.transpose()).trace()*
							(Yd.adjoint() * TYd)(i, j) - 2.0*(TYe.conjugate() * Ye.transpose()).trace()*
							(Yd.adjoint() * TYd)(i, j) +
							(8.0*std::pow(g(0), 2.0)*MHu2*(Yu.adjoint() * Yu)(i, j))/5.0
							- 12.0*MHu2*(Yu * Yu.adjoint()).trace()*
							(Yu.adjoint() * Yu)(i, j) - 6.0*(TYu.conjugate() * TYu.transpose()).trace()*
							(Yu.adjoint() * Yu)(i, j) - 6.0*(Mq2 * Yu.adjoint() * Yu).trace()*
							(Yu.adjoint() * Yu)(i, j) - 6.0*(Mu2 * Yu * Yu.adjoint()).trace()*
							(Yu.adjoint() * Yu)(i, j) +
							(std::pow(g(0), 2.0)*std::conj(M(0))*((597.0*std::pow(g(0), 2.0)*M(0)
							+ 80.0*std::pow(g(2), 2.0)*(2.0*M(0) + M(2)) +
							45.0*std::pow(g(1), 2.0)*(2.0*M(0) + M(1)))*KroneckerDelta(i, j) +
							180.0*(2.0*M(0)*(Yd.adjoint() * Yd)(i, j) - (Yd.adjoint() * TYd)(i, j)
							+ 4.0*M(0)*(Yu.adjoint() * Yu)(i, j) -
							2.0*(Yu.adjoint() * TYu)(i, j))))/225.0 -
							6.0*(TYu.conjugate() * Yu.transpose()).trace()*(Yu.adjoint() * TYu)(i, j) -
							(4.0*std::pow(g(0), 2.0)*M(0)*(TYd.adjoint() * Yd)(i, j))/5.0 -
							6.0*(Yd.adjoint() * TYd).trace()*(TYd.adjoint() * Yd)(i, j) -
							2.0*(Ye.adjoint() * TYe).trace()*(TYd.adjoint() * Yd)(i, j) +
							(4.0*std::pow(g(0), 2.0)*(TYd.adjoint() * TYd)(i, j))/5.0 -
							6.0*(Yd * Yd.adjoint()).trace()*(TYd.adjoint() * TYd)(i, j) -
							2.0*(Ye * Ye.adjoint()).trace()*(TYd.adjoint() * TYd)(i, j) -
							(8.0*std::pow(g(0), 2.0)*M(0)*(TYu.adjoint() * Yu)(i, j))/5.0 -
							6.0*(Yu.adjoint() * TYu).trace()*(TYu.adjoint() * Yu)(i, j) +
							(8.0*std::pow(g(0), 2.0)*(TYu.adjoint() * TYu)(i, j))/5.0 -
							6.0*(Yu * Yu.adjoint()).trace()*(TYu.adjoint() * TYu)(i, j) +
							(2.0*std::pow(g(0), 2.0)*(Mq2 * Yd.adjoint() * Yd)(i, j))/5.0 -
							3.0*(Yd * Yd.adjoint()).trace()*(Mq2 * Yd.adjoint() * Yd)(i, j) -
							(Ye * Ye.adjoint()).trace()*(Mq2 * Yd.adjoint() * Yd)(i, j) +
							(4.0*std::pow(g(0), 2.0)*(Mq2 * Yu.adjoint() * Yu)(i, j))/5.0 -
							3.0*(Yu * Yu.adjoint()).trace()*(Mq2 * Yu.adjoint() * Yu)(i, j) +
							(4.0*std::pow(g(0), 2.0)*(Yd.adjoint() * Md2 * Yd)(i, j))/5.0 -
							6.0*(Yd * Yd.adjoint()).trace()*(Yd.adjoint() * Md2 * Yd)(i, j) -
							2.0*(Ye * Ye.adjoint()).trace()*(Yd.adjoint() * Md2 * Yd)(i, j) +
							(2.0*std::pow(g(0), 2.0)*(Yd.adjoint() * Yd * Mq2)(i, j))/5.0 -
							3.0*(Yd * Yd.adjoint()).trace()*(Yd.adjoint() * Yd * Mq2)(i, j) -
							(Ye * Ye.adjoint()).trace()*(Yd.adjoint() * Yd * Mq2)(i, j) +
							(8.0*std::pow(g(0), 2.0)*(Yu.adjoint() * Mu2 * Yu)(i, j))/5.0 -
							6.0*(Yu * Yu.adjoint()).trace()*(Yu.adjoint() * Mu2 * Yu)(i, j) +
							(4.0*std::pow(g(0), 2.0)*(Yu.adjoint() * Yu * Mq2)(i, j))/5.0 -
							3.0*(Yu * Yu.adjoint()).trace()*(Yu.adjoint() * Yu * Mq2)(i, j) -
							8.0*MHd2*(Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j) -
							4.0*(Yd.adjoint() * Yd * TYd.adjoint() * TYd)(i, j) -
							4.0*(Yd.adjoint() * TYd * TYd.adjoint() * Yd)(i, j) -
							8.0*MHu2*(Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j) -
							4.0*(Yu.adjoint() * Yu * TYu.adjoint() * TYu)(i, j) -
							4.0*(Yu.adjoint() * TYu * TYu.adjoint() * Yu)(i, j) -
							4.0*(TYd.adjoint() * Yd * Yd.adjoint() * TYd)(i, j) -
							4.0*(TYd.adjoint() * TYd * Yd.adjoint() * Yd)(i, j) -
							4.0*(TYu.adjoint() * Yu * Yu.adjoint() * TYu)(i, j) -
							4.0*(TYu.adjoint() * TYu * Yu.adjoint() * Yu)(i, j) -
							2.0*(Mq2 * Yd.adjoint() * Yd * Yd.adjoint() * Yd)(i, j) -
							2.0*(Mq2 * Yu.adjoint() * Yu * Yu.adjoint() * Yu)(i, j) -
							4.0*(Yd.adjoint() * Md2 * Yd * Yd.adjoint() * Yd)(i, j) -
							4.0*(Yd.adjoint() * Yd * Mq2 * Yd.adjoint() * Yd)(i, j) -
							4.0*(Yd.adjoint() * Yd * Yd.adjoint() * Md2 * Yd)(i, j) -
							2.0*(Yd.adjoint() * Yd * Yd.adjoint() * Yd * Mq2)(i, j) -
							4.0*(Yu.adjoint() * Mu2 * Yu * Yu.adjoint() * Yu)(i, j) -
							4.0*(Yu.adjoint() * Yu * Mq2 * Yu.adjoint() * Yu)(i, j) -
							4.0*(Yu.adjoint() * Yu * Yu.adjoint() * Mu2 * Yu)(i, j) -
							2.0*(Yu.adjoint() * Yu * Yu.adjoint() * Yu * Mq2)(i, j);

			BetaMl2l1(i, j) = (-6.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0))*KroneckerDelta(i, j))/5.0 -
							6.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))*KroneckerDelta(i, j) -
							std::sqrt(3.0/5.0)*g(0)*KroneckerDelta(i, j)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2
							+ (Md2).trace() + (Me2).trace() -
							(Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace())) +
							2.0*MHd2*(Ye.adjoint() * Ye)(i, j) + 2.0*(TYe.adjoint() * TYe)(i, j) +
							(Ml2 * Ye.adjoint() * Ye)(i, j) + 2.0*(Ye.adjoint() * Me2 * Ye)(i, j) +
							(Ye.adjoint() * Ye * Ml2)(i, j);
			BetaMl2l2(i, j) = (3.0*std::pow(g(1), 2.0)*(55.0*std::pow(g(1), 2.0)*M(1)
							+ 3.0*std::pow(g(0), 2.0)*(M(0) + 2.0*M(1)))*std::conj(M(1))*
							KroneckerDelta(i, j))/5.0 + 6.0*std::pow(g(1), 4.0)*KroneckerDelta(i, j)*((MHd2 + MHu2
							+ (Ml2).trace() + 3.0*(Mq2).trace())/2.0) +
							(6.0*std::pow(g(0), 2.0)*KroneckerDelta(i, j)*((std::pow(g(0), 2.0)*(3.0*MHd2
							+ 3.0*MHu2 + 2.0*(Md2).trace() + 6.0*(Me2).trace() +
							3.0*(Ml2).trace() + (Mq2).trace() + 8.0*(Mu2).trace()))/10.0))/5.0 -
							4.0*std::sqrt(3.0/5.0)*g(0)*KroneckerDelta(i, j)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
							- 45.0*std::pow(g(1), 2.0)*MHd2 + 9.0*std::pow(g(0), 2.0)*MHu2
							+ 45.0*std::pow(g(1), 2.0)*MHu2 + 4.0*(std::pow(g(0), 2.0)
							+ 20.0*std::pow(g(2), 2.0))*(Md2).trace()
							+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
							9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
							+ std::pow(g(0), 2.0)*(Mq2).trace() +
							45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
							- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
							160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
							30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
							60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
							- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
							60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
							+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
							120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
							- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
							(20.0*std::sqrt(15.0))) +
							(12.0*std::pow(g(0), 2.0)*MHd2*(Ye.adjoint() * Ye)(i, j))/5.0 -
							12.0*MHd2*(Yd * Yd.adjoint()).trace()*(Ye.adjoint() * Ye)(i, j) -
							4.0*MHd2*(Ye * Ye.adjoint()).trace()*(Ye.adjoint() * Ye)(i, j) -
							6.0*(TYd.conjugate() * TYd.transpose()).trace()*(Ye.adjoint() * Ye)(i, j) -
							2.0*(TYe.conjugate() * TYe.transpose()).trace()*(Ye.adjoint() * Ye)(i, j) -
							6.0*(Md2 * Yd * Yd.adjoint()).trace()*(Ye.adjoint() * Ye)(i, j) -
							2.0*(Me2 * Ye * Ye.adjoint()).trace()*(Ye.adjoint() * Ye)(i, j) -
							2.0*(Ml2 * Ye.adjoint() * Ye).trace()*(Ye.adjoint() * Ye)(i, j) -
							6.0*(Mq2 * Yd.adjoint() * Yd).trace()*(Ye.adjoint() * Ye)(i, j) +
							(3.0*std::pow(g(0), 2.0)*std::conj(M(0))*(3.0*(69.0*std::pow(g(0), 2.0)*M(0)
							+ 5.0*std::pow(g(1), 2.0)*(2.0*M(0) + M(1)))*
							KroneckerDelta(i, j) + 40.0*M(0)*(Ye.adjoint() * Ye)(i, j) -
							20.0*(Ye.adjoint() * TYe)(i, j)))/25.0 - 6.0*(TYd.conjugate() * Yd.transpose()).trace()*
							(Ye.adjoint() * TYe)(i, j) - 2.0*(TYe.conjugate() * Ye.transpose()).trace()*
							(Ye.adjoint() * TYe)(i, j) -
							(12.0*std::pow(g(0), 2.0)*M(0)*(TYe.adjoint() * Ye)(i, j))/5.0 -
							6.0*(Yd.adjoint() * TYd).trace()*(TYe.adjoint() * Ye)(i, j) -
							2.0*(Ye.adjoint() * TYe).trace()*(TYe.adjoint() * Ye)(i, j) +
							(12.0*std::pow(g(0), 2.0)*(TYe.adjoint() * TYe)(i, j))/5.0 -
							6.0*(Yd * Yd.adjoint()).trace()*(TYe.adjoint() * TYe)(i, j) -
							2.0*(Ye * Ye.adjoint()).trace()*(TYe.adjoint() * TYe)(i, j) +
							(6.0*std::pow(g(0), 2.0)*(Ml2 * Ye.adjoint() * Ye)(i, j))/5.0 -
							3.0*(Yd * Yd.adjoint()).trace()*(Ml2 * Ye.adjoint() * Ye)(i, j) -
							(Ye * Ye.adjoint()).trace()*(Ml2 * Ye.adjoint() * Ye)(i, j) +
							(12.0*std::pow(g(0), 2.0)*(Ye.adjoint() * Me2 * Ye)(i, j))/5.0 -
							6.0*(Yd * Yd.adjoint()).trace()*(Ye.adjoint() * Me2 * Ye)(i, j) -
							2.0*(Ye * Ye.adjoint()).trace()*(Ye.adjoint() * Me2 * Ye)(i, j) +
							(6.0*std::pow(g(0), 2.0)*(Ye.adjoint() * Ye * Ml2)(i, j))/5.0 -
							3.0*(Yd * Yd.adjoint()).trace()*(Ye.adjoint() * Ye * Ml2)(i, j) -
							(Ye * Ye.adjoint()).trace()*(Ye.adjoint() * Ye * Ml2)(i, j) -
							8.0*MHd2*(Ye.adjoint() * Ye * Ye.adjoint() * Ye)(i, j) -
							4.0*(Ye.adjoint() * Ye * TYe.adjoint() * TYe)(i, j) -
							4.0*(Ye.adjoint() * TYe * TYe.adjoint() * Ye)(i, j) -
							4.0*(TYe.adjoint() * Ye * Ye.adjoint() * TYe)(i, j) -
							4.0*(TYe.adjoint() * TYe * Ye.adjoint() * Ye)(i, j) -
							2.0*(Ml2 * Ye.adjoint() * Ye * Ye.adjoint() * Ye)(i, j) -
							4.0*(Ye.adjoint() * Me2 * Ye * Ye.adjoint() * Ye)(i, j) -
							4.0*(Ye.adjoint() * Ye * Ml2 * Ye.adjoint() * Ye)(i, j) -
							4.0*(Ye.adjoint() * Ye * Ye.adjoint() * Me2 * Ye)(i, j) -
							2.0*(Ye.adjoint() * Ye * Ye.adjoint() * Ye * Ml2)(i, j);

			BetaMd2l1(i, j) = (-8.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0))*KroneckerDelta(i, j))/15.0 -
							(32.0*std::pow(g(2), 2.0)*M(2)*std::conj(M(2))*KroneckerDelta(i, j))/3.0 +
							(2.0*g(0)*KroneckerDelta(i, j)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace()
							+ (Me2).trace() - (Ml2).trace() + (Mq2).trace()
							- 2.0*(Mu2).trace())))/std::sqrt(15.0) +
							4.0*MHd2*(Yd * Yd.adjoint())(i, j) + 4.0*(TYd * TYd.adjoint())(i, j) +
							2.0*(Md2 * Yd * Yd.adjoint())(i, j) + 4.0*(Yd * Mq2 * Yd.adjoint())(i, j) +
							2.0*(Yd * Yd.adjoint() * Md2)(i, j);
			BetaMd2l2(i, j) = (64.0*std::pow(g(2), 2.0)*(-30.0*std::pow(g(2), 2.0)*M(2)
							+ std::pow(g(0), 2.0)*(M(0) + 2.0*M(2)))*std::conj(M(2))*
							KroneckerDelta(i, j))/45.0 + (32.0*std::pow(g(2), 4.0)*KroneckerDelta(i, j)*(((Md2).trace()
							+ 2.0*(Mq2).trace() + (Mu2).trace())/2.0))/3.0 +
							(8.0*std::pow(g(0), 2.0)*KroneckerDelta(i, j)*((std::pow(g(0), 2.0)*(3.0*MHd2
							+ 3.0*MHu2 + 2.0*(Md2).trace() + 6.0*(Me2).trace() +
							3.0*(Ml2).trace() + (Mq2).trace() + 8.0*(Mu2).trace()))/10.0))/15.0 +
							(8.0*g(0)*KroneckerDelta(i, j)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
							- 45.0*std::pow(g(1), 2.0)*MHd2
							+ 9.0*std::pow(g(0), 2.0)*MHu2 + 45.0*std::pow(g(1), 2.0)*MHu2 +
							4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Md2).trace()
							+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
							9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
							+ std::pow(g(0), 2.0)*(Mq2).trace() +
							45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
							- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
							160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
							30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
							60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
							- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
							60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
							+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
							120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
							- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
							(20.0*std::sqrt(15.0))))/std::sqrt(15.0) +
							(4.0*std::pow(g(0), 2.0)*MHd2*(Yd * Yd.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*MHd2*(Yd * Yd.adjoint())(i, j)
							+ 24.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))*
							(Yd * Yd.adjoint())(i, j) - 24.0*MHd2*(Yd * Yd.adjoint()).trace()*
							(Yd * Yd.adjoint())(i, j) - 8.0*MHd2*(Ye * Ye.adjoint()).trace()*
							(Yd * Yd.adjoint())(i, j) - 12.0*(TYd.conjugate() * TYd.transpose()).trace()*
							(Yd * Yd.adjoint())(i, j) - 4.0*(TYe.conjugate() * TYe.transpose()).trace()*
							(Yd * Yd.adjoint())(i, j) - 12.0*(Md2 * Yd * Yd.adjoint()).trace()*
							(Yd * Yd.adjoint())(i, j) - 4.0*(Me2 * Ye * Ye.adjoint()).trace()*
							(Yd * Yd.adjoint())(i, j) - 4.0*(Ml2 * Ye.adjoint() * Ye).trace()*
							(Yd * Yd.adjoint())(i, j) - 12.0*(Mq2 * Yd.adjoint() * Yd).trace()*
							(Yd * Yd.adjoint())(i, j) -
							(4.0*std::pow(g(0), 2.0)*M(0)*(Yd * TYd.adjoint())(i, j))/5.0 -
							12.0*std::pow(g(1), 2.0)*M(1)*(Yd * TYd.adjoint())(i, j) -
							12.0*(Yd.adjoint() * TYd).trace()*(Yd * TYd.adjoint())(i, j) -
							4.0*(Ye.adjoint() * TYe).trace()*(Yd * TYd.adjoint())(i, j) +
							(4.0*std::pow(g(0), 2.0)*std::conj(M(0))*(2.0*(303.0*std::pow(g(0), 2.0)*M(0)
							+ 40.0*std::pow(g(2), 2.0)*(2.0*M(0) + M(2)))*
							KroneckerDelta(i, j) + 90.0*M(0)*(Yd * Yd.adjoint())(i, j) -
							45.0*(TYd * Yd.adjoint())(i, j)))/225.0 -
							12.0*std::pow(g(1), 2.0)*std::conj(M(1))*(TYd * Yd.adjoint())(i, j) -
							12.0*(TYd.conjugate() * Yd.transpose()).trace()*(TYd * Yd.adjoint())(i, j) -
							4.0*(TYe.conjugate() * Ye.transpose()).trace()*(TYd * Yd.adjoint())(i, j) +
							(4.0*std::pow(g(0), 2.0)*(TYd * TYd.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*(TYd * TYd.adjoint())(i, j) - 12.0*(Yd * Yd.adjoint()).trace()*
							(TYd * TYd.adjoint())(i, j) - 4.0*(Ye * Ye.adjoint()).trace()*
							(TYd * TYd.adjoint())(i, j) +
							(2.0*std::pow(g(0), 2.0)*(Md2 * Yd * Yd.adjoint())(i, j))/5.0 +
							6.0*std::pow(g(1), 2.0)*(Md2 * Yd * Yd.adjoint())(i, j) - 6.0*(Yd * Yd.adjoint()).trace()*
							(Md2 * Yd * Yd.adjoint())(i, j) - 2.0*(Ye * Ye.adjoint()).trace()*
							(Md2 * Yd * Yd.adjoint())(i, j) +
							(4.0*std::pow(g(0), 2.0)*(Yd * Mq2 * Yd.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*(Yd * Mq2 * Yd.adjoint())(i, j) - 12.0*(Yd * Yd.adjoint()).trace()*
							(Yd * Mq2 * Yd.adjoint())(i, j) - 4.0*(Ye * Ye.adjoint()).trace()*
							(Yd * Mq2 * Yd.adjoint())(i, j) +
							(2.0*std::pow(g(0), 2.0)*(Yd * Yd.adjoint() * Md2)(i, j))/5.0 +
							6.0*std::pow(g(1), 2.0)*(Yd * Yd.adjoint() * Md2)(i, j) - 6.0*(Yd * Yd.adjoint()).trace()*
							(Yd * Yd.adjoint() * Md2)(i, j) - 2.0*(Ye * Ye.adjoint()).trace()*
							(Yd * Yd.adjoint() * Md2)(i, j) -
							8.0*MHd2*(Yd * Yd.adjoint() * Yd * Yd.adjoint())(i, j) -
							4.0*(Yd * Yd.adjoint() * TYd * TYd.adjoint())(i, j) -
							4.0*MHd2*(Yd * Yu.adjoint() * Yu * Yd.adjoint())(i, j) -
							4.0*MHu2*(Yd * Yu.adjoint() * Yu * Yd.adjoint())(i, j) -
							4.0*(Yd * Yu.adjoint() * TYu * TYd.adjoint())(i, j) -
							4.0*(Yd * TYd.adjoint() * TYd * Yd.adjoint())(i, j) -
							4.0*(Yd * TYu.adjoint() * TYu * Yd.adjoint())(i, j) -
							4.0*(TYd * Yd.adjoint() * Yd * TYd.adjoint())(i, j) -
							4.0*(TYd * Yu.adjoint() * Yu * TYd.adjoint())(i, j) -
							4.0*(TYd * TYd.adjoint() * Yd * Yd.adjoint())(i, j) -
							4.0*(TYd * TYu.adjoint() * Yu * Yd.adjoint())(i, j) -
							2.0*(Md2 * Yd * Yd.adjoint() * Yd * Yd.adjoint())(i, j) -
							2.0*(Md2 * Yd * Yu.adjoint() * Yu * Yd.adjoint())(i, j) -
							4.0*(Yd * Mq2 * Yd.adjoint() * Yd * Yd.adjoint())(i, j) -
							4.0*(Yd * Mq2 * Yu.adjoint() * Yu * Yd.adjoint())(i, j) -
							4.0*(Yd * Yd.adjoint() * Md2 * Yd * Yd.adjoint())(i, j) -
							4.0*(Yd * Yd.adjoint() * Yd * Mq2 * Yd.adjoint())(i, j) -
							2.0*(Yd * Yd.adjoint() * Yd * Yd.adjoint() * Md2)(i, j) -
							4.0*(Yd * Yu.adjoint() * Mu2 * Yu * Yd.adjoint())(i, j) -
							4.0*(Yd * Yu.adjoint() * Yu * Mq2 * Yd.adjoint())(i, j) -
							2.0*(Yd * Yu.adjoint() * Yu * Yd.adjoint() * Md2)(i, j);

			BetaMu2l1(i, j) = (-32.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0))*KroneckerDelta(i, j))/15.0 -
							(32.0*std::pow(g(2), 2.0)*M(2)*std::conj(M(2))*KroneckerDelta(i, j))/3.0 -
							(4.0*g(0)*KroneckerDelta(i, j)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace()
							+ (Me2).trace() - (Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace())))/std::sqrt(15.0) +
							4.0*MHu2*(Yu * Yu.adjoint())(i, j) + 4.0*(TYu * TYu.adjoint())(i, j) +
							2.0*(Mu2 * Yu * Yu.adjoint())(i, j) + 4.0*(Yu * Mq2 * Yu.adjoint())(i, j) +
							2.0*(Yu * Yu.adjoint() * Mu2)(i, j);
			BetaMu2l2(i, j) = (-128.0*std::pow(g(2), 2.0)*(15.0*std::pow(g(2), 2.0)*M(2)
							- 2.0*std::pow(g(0), 2.0)*(M(0) + 2.0*M(2)))*std::conj(M(2))*
							KroneckerDelta(i, j))/45.0 + (32.0*std::pow(g(2), 4.0)*KroneckerDelta(i, j)*(((Md2).trace()
							+ 2.0*(Mq2).trace() + (Mu2).trace())/2.0))/3.0 +
							(32.0*std::pow(g(0), 2.0)*KroneckerDelta(i, j)*((std::pow(g(0), 2.0)*(3.0*MHd2 + 3.0*MHu2
							+ 2.0*(Md2).trace() + 6.0*(Me2).trace() + 3.0*(Ml2).trace() + (Mq2).trace()
							+ 8.0*(Mu2).trace()))/10.0))/15.0
							- (16.0*g(0)*KroneckerDelta(i, j)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
							- 45.0*std::pow(g(1), 2.0)*MHd2 + 9.0*std::pow(g(0), 2.0)*MHu2
							+ 45.0*std::pow(g(1), 2.0)*MHu2 + 4.0*(std::pow(g(0), 2.0)
							+ 20.0*std::pow(g(2), 2.0))*(Md2).trace()
							+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() -
							9.0*std::pow(g(0), 2.0)*(Ml2).trace() - 45.0*std::pow(g(1), 2.0)*(Ml2).trace()
							+ std::pow(g(0), 2.0)*(Mq2).trace() +
							45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
							- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() -
							160.0*std::pow(g(2), 2.0)*(Mu2).trace() + 90.0*MHd2*(Yd * (Yd).adjoint()).trace() +
							30.0*MHd2*(Ye * (Ye).adjoint()).trace() - 90.0*MHu2*(Yu * (Yu).adjoint()).trace() -
							60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
							- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace() -
							60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
							+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace() +
							120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
							- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/
							(20.0*std::sqrt(15.0))))/std::sqrt(15.0) -
							(4.0*std::pow(g(0), 2.0)*MHu2*(Yu * Yu.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*MHu2*(Yu * Yu.adjoint())(i, j)
							+ 24.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))*
							(Yu * Yu.adjoint())(i, j) - 24.0*MHu2*(Yu * Yu.adjoint()).trace()*
							(Yu * Yu.adjoint())(i, j) - 12.0*(TYu.conjugate() * TYu.transpose()).trace()*
							(Yu * Yu.adjoint())(i, j) - 12.0*(Mq2 * Yu.adjoint() * Yu).trace()*
							(Yu * Yu.adjoint())(i, j) - 12.0*(Mu2 * Yu * Yu.adjoint()).trace()*
							(Yu * Yu.adjoint())(i, j) +
							(4.0*std::pow(g(0), 2.0)*M(0)*(Yu * TYu.adjoint())(i, j))/5.0 -
							12.0*std::pow(g(1), 2.0)*M(1)*(Yu * TYu.adjoint())(i, j) -
							12.0*(Yu.adjoint() * TYu).trace()*(Yu * TYu.adjoint())(i, j) -
							12.0*std::pow(g(1), 2.0)*std::conj(M(1))*(TYu * Yu.adjoint())(i, j) -
							12.0*(TYu.conjugate() * Yu.transpose()).trace()*(TYu * Yu.adjoint())(i, j) +
							(4.0*std::pow(g(0), 2.0)*std::conj(M(0))*(8.0*(321.0*std::pow(g(0), 2.0)*M(0)
							+ 40.0*std::pow(g(2), 2.0)*(2.0*M(0) + M(2)))*
							KroneckerDelta(i, j) + 45.0*(-2.0*M(0)*(Yu * Yu.adjoint())(i, j) +
							(TYu * Yu.adjoint())(i, j))))/225.0 -
							(4.0*std::pow(g(0), 2.0)*(TYu * TYu.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*(TYu * TYu.adjoint())(i, j) - 12.0*(Yu * Yu.adjoint()).trace()*
							(TYu * TYu.adjoint())(i, j) -
							(2.0*std::pow(g(0), 2.0)*(Mu2 * Yu * Yu.adjoint())(i, j))/5.0 +
							6.0*std::pow(g(1), 2.0)*(Mu2 * Yu * Yu.adjoint())(i, j) - 6.0*(Yu * Yu.adjoint()).trace()*
							(Mu2 * Yu * Yu.adjoint())(i, j) -
							(4.0*std::pow(g(0), 2.0)*(Yu * Mq2 * Yu.adjoint())(i, j))/5.0 +
							12.0*std::pow(g(1), 2.0)*(Yu * Mq2 * Yu.adjoint())(i, j) - 12.0*(Yu * Yu.adjoint()).trace()*
							(Yu * Mq2 * Yu.adjoint())(i, j) -
							(2.0*std::pow(g(0), 2.0)*(Yu * Yu.adjoint() * Mu2)(i, j))/5.0 +
							6.0*std::pow(g(1), 2.0)*(Yu * Yu.adjoint() * Mu2)(i, j) - 6.0*(Yu * Yu.adjoint()).trace()*
							(Yu * Yu.adjoint() * Mu2)(i, j) -
							4.0*MHd2*(Yu * Yd.adjoint() * Yd * Yu.adjoint())(i, j) -
							4.0*MHu2*(Yu * Yd.adjoint() * Yd * Yu.adjoint())(i, j) -
							4.0*(Yu * Yd.adjoint() * TYd * TYu.adjoint())(i, j) -
							8.0*MHu2*(Yu * Yu.adjoint() * Yu * Yu.adjoint())(i, j) -
							4.0*(Yu * Yu.adjoint() * TYu * TYu.adjoint())(i, j) -
							4.0*(Yu * TYd.adjoint() * TYd * Yu.adjoint())(i, j) -
							4.0*(Yu * TYu.adjoint() * TYu * Yu.adjoint())(i, j) -
							4.0*(TYu * Yd.adjoint() * Yd * TYu.adjoint())(i, j) -
							4.0*(TYu * Yu.adjoint() * Yu * TYu.adjoint())(i, j) -
							4.0*(TYu * TYd.adjoint() * Yd * Yu.adjoint())(i, j) -
							4.0*(TYu * TYu.adjoint() * Yu * Yu.adjoint())(i, j) -
							2.0*(Mu2 * Yu * Yd.adjoint() * Yd * Yu.adjoint())(i, j) -
							2.0*(Mu2 * Yu * Yu.adjoint() * Yu * Yu.adjoint())(i, j) -
							4.0*(Yu * Mq2 * Yd.adjoint() * Yd * Yu.adjoint())(i, j) -
							4.0*(Yu * Mq2 * Yu.adjoint() * Yu * Yu.adjoint())(i, j) -
							4.0*(Yu * Yd.adjoint() * Md2 * Yd * Yu.adjoint())(i, j) -
							4.0*(Yu * Yd.adjoint() * Yd * Mq2 * Yu.adjoint())(i, j) -
							2.0*(Yu * Yd.adjoint() * Yd * Yu.adjoint() * Mu2)(i, j) -
							4.0*(Yu * Yu.adjoint() * Mu2 * Yu * Yu.adjoint())(i, j) -
							4.0*(Yu * Yu.adjoint() * Yu * Mq2 * Yu.adjoint())(i, j) -
							2.0*(Yu * Yu.adjoint() * Yu * Yu.adjoint() * Mu2)(i, j);

			BetaMe2l1(i, j) = (-24.0*std::pow(g(0), 2.0)*M(0)*std::conj(M(0))*KroneckerDelta(i, j))/5.0 +
							2.0*std::sqrt(3.0/5.0)*g(0)*KroneckerDelta(i, j)*(std::sqrt(3.0/5.0)*g(0)*(-MHd2 + MHu2 + (Md2).trace()
							+ (Me2).trace() - (Ml2).trace() + (Mq2).trace() - 2.0*(Mu2).trace()))
							+ 2.0*(2.0*MHd2*(Ye * Ye.adjoint())(i, j)
							+ 2.0*(TYe * TYe.adjoint())(i, j) + (Me2 * Ye * Ye.adjoint())(i, j)
							+ 2.0*(Ye * Ml2 * Ye.adjoint())(i, j) + (Ye * Ye.adjoint() * Me2)(i, j));
			BetaMe2l2(i, j) = (2.0*(20.0*g(0)*KroneckerDelta(i, j)*(3.0*g(0)*((std::pow(g(0), 2.0)*(3.0*MHd2
                            + 3.0*MHu2 + 2.0*(Md2).trace() + 6.0*(Me2).trace() + 3.0*(Ml2).trace()
							+ (Mq2).trace() + 8.0*(Mu2).trace()))/10.0)
							+ std::sqrt(15.0)*((g(0)*(-9.0*std::pow(g(0), 2.0)*MHd2
							- 45.0*std::pow(g(1), 2.0)*MHd2 + 9.0*std::pow(g(0), 2.0)*MHu2
							+ 45.0*std::pow(g(1), 2.0)*MHu2
							+ 4.0*(std::pow(g(0), 2.0) + 20.0*std::pow(g(2), 2.0))*(Md2).trace()
							+ 36.0*std::pow(g(0), 2.0)*(Me2).trace() - 9.0*std::pow(g(0), 2.0)*(Ml2).trace()
							- 45.0*std::pow(g(1), 2.0)*(Ml2).trace() + std::pow(g(0), 2.0)*(Mq2).trace()
							+ 45.0*std::pow(g(1), 2.0)*(Mq2).trace() + 80.0*std::pow(g(2), 2.0)*(Mq2).trace()
							- 32.0*std::pow(g(0), 2.0)*(Mu2).trace() - 160.0*std::pow(g(2), 2.0)*(Mu2).trace()
							+ 90.0*MHd2*(Yd * (Yd).adjoint()).trace() + 30.0*MHd2*(Ye * (Ye).adjoint()).trace()
							- 90.0*MHu2*(Yu * (Yu).adjoint()).trace() - 60.0*(Yd * (Yd).adjoint() * (Md2).conjugate()).trace()
							- 30.0*(Yd * (Mq2).conjugate() * (Yd).adjoint()).trace()
							- 60.0*(Ye * (Ye).adjoint() * (Me2).conjugate()).trace()
							+ 30.0*(Ye * (Ml2).conjugate() * (Ye).adjoint()).trace()
							+ 120.0*(Yu * (Yu).adjoint() * (Mu2).conjugate()).trace()
							- 30.0*(Yu * (Mq2).conjugate() * (Yu).adjoint()).trace()))/(20.0*std::sqrt(15.0))))
							+ 6.0*std::conj(M(0))*(234.0*std::pow(g(0), 4.0)*M(0)*KroneckerDelta(i, j)
							+ 5.0*std::pow(g(0), 2.0)*(-2.0*M(0)*(Ye * Ye.adjoint())(i, j)
							+ (TYe * Ye.adjoint())(i, j)))
							- 5.0*(2.0*(3.0*std::pow(g(0), 2.0)*MHd2 - 15.0*std::pow(g(1), 2.0)*MHd2
							- 30.0*std::pow(g(1), 2.0)*M(1)*std::conj(M(1))
							+ 30.0*MHd2*(Yd * Yd.adjoint()).trace() + 10.0*MHd2*(Ye * Ye.adjoint()).trace()
							+ 15.0*(TYd.conjugate() * TYd.transpose()).trace()
							+ 5.0*(TYe.conjugate() * TYe.transpose()).trace() + 15.0*(Md2 * Yd * Yd.adjoint()).trace()
							+ 5.0*(Me2 * Ye * Ye.adjoint()).trace() + 5.0*(Ml2 * Ye.adjoint() * Ye).trace()
							+ 15.0*(Mq2 * Yd.adjoint() * Yd).trace())*(Ye * Ye.adjoint())(i, j)
							+ (-6.0*std::pow(g(0), 2.0)*M(0) + 30.0*std::pow(g(1), 2.0)*M(1)
							+ 30.0*(Yd.adjoint() * TYd).trace()
							+ 10.0*(Ye.adjoint() * TYe).trace())*(Ye * TYe.adjoint())(i, j)
							+ 30.0*std::pow(g(1), 2.0)*std::conj(M(1))*(TYe * Ye.adjoint())(i, j)
							+ 30.0*(TYd.conjugate() * Yd.transpose()).trace()*(TYe * Ye.adjoint())(i, j)
							+ 10.0*(TYe.conjugate() * Ye.transpose()).trace()*(TYe * Ye.adjoint())(i, j) +
							6.0*std::pow(g(0), 2.0)*(TYe * TYe.adjoint())(i, j) -
							30.0*std::pow(g(1), 2.0)*(TYe * TYe.adjoint())(i, j)
							+ 30.0*(Yd * Yd.adjoint()).trace()*
							(TYe * TYe.adjoint())(i, j) + 10.0*(Ye * Ye.adjoint()).trace()*
							(TYe * TYe.adjoint())(i, j)
							+ 3.0*std::pow(g(0), 2.0)*(Me2 * Ye * Ye.adjoint())(i, j)
							- 15.0*std::pow(g(1), 2.0)*(Me2 * Ye * Ye.adjoint())(i, j) +
							15.0*(Yd * Yd.adjoint()).trace()*(Me2 * Ye * Ye.adjoint())(i, j) +
							5.0*(Ye * Ye.adjoint()).trace()*(Me2 * Ye * Ye.adjoint())(i, j) +
							6.0*std::pow(g(0), 2.0)*(Ye * Ml2 * Ye.adjoint())(i, j) -
							30.0*std::pow(g(1), 2.0)*(Ye * Ml2 * Ye.adjoint())(i, j)
							+ 30.0*(Yd * Yd.adjoint()).trace()*
							(Ye * Ml2 * Ye.adjoint())(i, j) + 10.0*(Ye * Ye.adjoint()).trace()*
							(Ye * Ml2 * Ye.adjoint())(i, j)
							+ 3.0*std::pow(g(0), 2.0)*(Ye * Ye.adjoint() * Me2)(i, j)
							- 15.0*std::pow(g(1), 2.0)*(Ye * Ye.adjoint() * Me2)(i, j) +
							15.0*(Yd * Yd.adjoint()).trace()*(Ye * Ye.adjoint() * Me2)(i, j) +
							5.0*(Ye * Ye.adjoint()).trace()*(Ye * Ye.adjoint() * Me2)(i, j) +
							20.0*MHd2*(Ye * Ye.adjoint() * Ye * Ye.adjoint())(i, j) +
							10.0*(Ye * Ye.adjoint() * TYe * TYe.adjoint())(i, j) +
							10.0*(Ye * TYe.adjoint() * TYe * Ye.adjoint())(i, j) +
							10.0*(TYe * Ye.adjoint() * Ye * TYe.adjoint())(i, j) +
							10.0*(TYe * TYe.adjoint() * Ye * Ye.adjoint())(i, j) +
							5.0*(Me2 * Ye * Ye.adjoint() * Ye * Ye.adjoint())(i, j) +
							10.0*(Ye * Ml2 * Ye.adjoint() * Ye * Ye.adjoint())(i, j) +
							10.0*(Ye * Ye.adjoint() * Me2 * Ye * Ye.adjoint())(i, j) +
							10.0*(Ye * Ye.adjoint() * Ye * Ml2 * Ye.adjoint())(i, j) +
							5.0*(Ye * Ye.adjoint() * Ye * Ye.adjoint() * Me2)(i, j))))/25.0;
		}
	}

	//Adding the 2 loop orders of the Beta functions:
	dg1_dt = LoopFactor * Betag1l1 + std::pow(LoopFactor, 2) * Betag1l2;
	dg2_dt = LoopFactor * Betag2l1 + std::pow(LoopFactor, 2) * Betag2l2;
	dg3_dt = LoopFactor * Betag3l1 + std::pow(LoopFactor, 2) * Betag3l2;

	dM1_dt = LoopFactor * BetaM1l1 + std::pow(LoopFactor, 2) * BetaM1l2;
	dM2_dt = LoopFactor * BetaM2l1 + std::pow(LoopFactor, 2) * BetaM2l2;
	dM3_dt = LoopFactor * BetaM3l1 + std::pow(LoopFactor, 2) * BetaM3l2;

	dMHu2_dt = LoopFactor * BetaMHu2l1 + std::pow(LoopFactor, 2) * BetaMHu2l2;
	dMHd2_dt = LoopFactor * BetaMHd2l1 + std::pow(LoopFactor, 2) * BetaMHd2l2;

	//dBmu_dt = LoopFactor * BetaBmul1 + std::pow(LoopFactor, 2) * BetaBmul2;
	dYd_dt = LoopFactor * BetaYdl1 + std::pow(LoopFactor, 2) * BetaYdl2;
	dYu_dt = LoopFactor * BetaYul1 + std::pow(LoopFactor, 2) * BetaYul2;
	dYe_dt = LoopFactor * BetaYel1 + std::pow(LoopFactor, 2) * BetaYel2;
	dTYd_dt = LoopFactor * BetaTYdl1 + std::pow(LoopFactor, 2) * BetaTYdl2;
	dTYu_dt = LoopFactor * BetaTYul1 + std::pow(LoopFactor, 2) * BetaTYul2;
	dTYe_dt = LoopFactor * BetaTYel1 + std::pow(LoopFactor, 2) * BetaTYel2;
	dMq2_dt = LoopFactor * BetaMq2l1 + std::pow(LoopFactor, 2) * BetaMq2l2;
	dMl2_dt = LoopFactor * BetaMl2l1 + std::pow(LoopFactor, 2) * BetaMl2l2;
	dMu2_dt = LoopFactor * BetaMu2l1 + std::pow(LoopFactor, 2) * BetaMu2l2;
	dMd2_dt = LoopFactor * BetaMd2l1 + std::pow(LoopFactor, 2) * BetaMd2l2;
	dMe2_dt = LoopFactor * BetaMe2l1 + std::pow(LoopFactor, 2) * BetaMe2l2;

	//storing the derivatives dy(t)/dt into one vector to facilitate integration:
	VectorXcd dydt(y.size());
	k = 0;
	for(int i=0; i<=2; i++) {
		for (int j=0; j<=2; j++)
		{
			dydt(k) = dYu_dt(i,j);
			dydt(9+k)  = dYd_dt(i,j);
			dydt(18+k) = dYe_dt(i,j);
			dydt(27+k) = dTYu_dt(i,j);
			dydt(36+k) = dTYd_dt(i,j);
			dydt(45+k) = dTYe_dt(i,j);
			dydt(54+k) = dMq2_dt(i,j);
			dydt(63+k) = dMl2_dt(i,j);
			dydt(72+k) = dMu2_dt(i,j);
			dydt(81+k) = dMd2_dt(i,j);
			dydt(90+k) = dMe2_dt(i,j);
			k++;
		}
	}

	dydt(99) = dM1_dt;
	dydt(100) = dM2_dt;
	dydt(101) = dM3_dt;
	dydt(102) = dg1_dt;
	dydt(103) = dg2_dt;
	dydt(104) = dg3_dt;
	dydt(105) = dMHu2_dt;
	dydt(106) = dMHd2_dt;
	//dydt(107) = dBmu_dt;

	return dydt;
}











































