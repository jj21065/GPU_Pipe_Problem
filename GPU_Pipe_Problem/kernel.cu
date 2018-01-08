
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
float myabs(float);
float mypow(float, int);
int main()
{
	// pipe parameter 
	float r01 = 1;
	float r16 = 5;
	float r12 = 1;
	float r03 = 5;
	float r34 = 3;

	float r24 = 1;
	float r25 = 1;
	float r56 = 3;
	float r3e = 1;
	float r47 = 1;

	float r57 = 1;
	float r78 = 2;
	float r89 = 3;
	float r69 = 3;

	float r8a = 1;
	float reb = 2;
	float rdg = 5;
	float ref = 5;

	float rab = 3;
	float rbc = 2;
	float r9g = 1;
	float rcf = 2;

	float rad = 2;
	float rdc = 1;
	float rgh = 1;
	float rfh = 2;

	float rs0 = 1;
	float r9EA = 1;
	float rhEB = 1;


	float Qs0 = Qs0 = sqrt(10.0 / rs0);
	float Q01 = Qs0;
	float Q16 = Qs0;
	float Q12 = 0;
	float Q03 = 0;
	float Q34 = 0;
	float Q24 = 0;
	float Q25 = 0;
	float Q56 = 0;
	float Q3e = 0;
	float Q47 = 0;

	float Q57 = 0;
	float Q78 = 0;
	float Q89 = 0;
	float Q69 = Qs0;

	float Q8a = 0;
	float Qeb = 0;
	float Qdg = 0;
	float Qef = 0;

	float Qab = 0;
	float Qbc = 0;
	float Q9g = 0;
	float Qcf = 0;

	float Qad = 0;
	float Qdc = 0;
	float Qgh = 0;
	float Qfh = 0;


	float Q9EA = Qs0;
	float QhEB = 0;


	float R[11] = { 0 };

	// initial conditions 


	int n = 2;
	int iter_no = 5000;
	int i;
	clock_t t1 = clock();
	printf("inital Qs0 = %g\n", Qs0);
	for (i = 0; i < iter_no; i++)
	{
		/// calculate adjust flow Q
		if (2 * (r01*abs(Q01) + r12*abs(Q12) + r24*abs(Q24) + r34*abs(Q34) + r03*abs(Q03)) > 0)
			R[0] = -(r01*Q01*abs(Q01) + r12*Q12*abs(Q12) + r24*Q24*abs(Q24) - r34*Q34*abs(Q34) - r03*Q03*abs(Q03)) / (2 * (r01*abs(Q01) + r12*abs(Q12) + r24*abs(Q24) + r34*abs(Q34) + r03*abs(Q03)));
		
		if (2 * (r16*abs(Q16) + r56*abs(Q56) + r25*abs(Q25) + r12*abs(Q12)) > 0)
			R[1] = -(r16*Q16*abs(Q16) - r56*Q56*abs(Q56) - r25*Q25*abs(Q25) - r12*Q12*abs(Q12)) / (2 * (r16*abs(Q16) + r56*abs(Q56) + r25*abs(Q25) + r12*abs(Q12)));
		
		if (2 * (r34*abs(Q34) + r47*abs(Q47) + r78*abs(Q78) + r8a*abs(Q8a) + rab*abs(Qab) + reb*abs(Qeb) + r3e*abs(Q3e)) > 0)
			R[2] = -(r34*Q34*abs(Q34) + r47*Q47*abs(Q47) + r78*Q78*abs(Q78) + r8a*Q8a*abs(Q8a) + rab*Qab*abs(Qab) - reb*Qeb*abs(Qeb) - r3e*Q3e*abs(Q3e)) / (2 * (r34*abs(Q34) + r47*abs(Q47) + r78*abs(Q78) + r8a*abs(Q8a) + rab*abs(Qab) + reb*abs(Qeb) + r3e*abs(Q3e)));
		
		if (2 * (r56*abs(Q56) + r69*abs(Q69) + r89*abs(Q89) + r78*abs(Q78) + r57*abs(Q57)) > 0)
			R[3] = -(r56*Q56*abs(Q56) + r69*Q69*abs(Q69) - r89*Q89*abs(Q89) - r78*Q78*abs(Q78) - r57*Q57*abs(Q57)) / (2 * (r56*abs(Q56) + r69*abs(Q69) + r89*abs(Q89) + r78*abs(Q78) + r57*abs(Q57)));
		
		if (2 * (r89*abs(Q89) + r9g*abs(Q9g) + rdg*abs(Qdg) + rad*abs(Qad) + r8a*abs(Q8a)) > 0)
			R[4] = -(r89*Q89*abs(Q89) - r9g*Q9g*abs(Q9g) - rdg*Qdg*abs(Qdg) - rad*Qad*abs(Qad) - r8a*Q8a*abs(Q8a)) / (2 * (r89*abs(Q89) + r9g*abs(Q9g) + rdg*abs(Qdg) + rad*abs(Qad) + r8a*abs(Q8a)));
		
		if (2 * (rad*abs(Qad) + rdc*abs(Qdc) + rbc*abs(Qbc) + rab*abs(Qab)) > 0)
			R[5] = -(rad*Qad*abs(Qad) + rdc*Qdc*abs(Qdc) - rbc*Qbc*abs(Qbc) - rab*Qab*abs(Qab)) / (2 * (rad*abs(Qad) + rdc*abs(Qdc) + rbc*abs(Qbc) + rab*abs(Qab)));
		
		if (2 * (reb*abs(Qeb) + rbc*abs(Qbc) + rcf*abs(Qcf) + ref*abs(Qef)) > 0)
			R[6] = -(reb*Qeb*abs(Qeb) + rbc*Qbc*abs(Qbc) + rcf*Qcf*abs(Qcf) - ref*Qef*abs(Qef)) / (2 * (reb*abs(Qeb) + rbc*abs(Qbc) + rcf*abs(Qcf) + ref*abs(Qef)));
		
		if (2 * (rdg*abs(Qdg) + rgh*abs(Qgh) + rfh*abs(Qfh) + rcf*abs(Qcf) + rdc*abs(Qdc)) > 0)
			R[7] = -(rdg*Qdg*abs(Qdg) + rgh*Qgh*abs(Qgh) - rfh*Qfh*abs(Qfh) - rcf*Qcf*abs(Qcf) - rdc*Qdc*abs(Qdc)) / (2 * (rdg*abs(Qdg) + rgh*abs(Qgh) + rfh*abs(Qfh) + rcf*abs(Qcf) + rdc*abs(Qdc)));
		
		if (2 * (r25*abs(Q25) + r57*abs(Q57) + r47*abs(Q47) + r24*abs(Q24)) > 0)
			R[8] = -(r25*Q25*abs(Q25) + r57*Q57*abs(Q57) - r47*Q47*abs(Q47) - r24*Q24*abs(Q24)) / (2 * (r25*abs(Q25) + r57*abs(Q57) + r47*abs(Q47) + r24*abs(Q24)));
		
		if (2 * (rs0*abs(Qs0) + r01*abs(Q01) + r16*abs(Q16) + r69*abs(Q69) + r9EA*abs(Q9EA)) > 0)
			R[9] = -(rs0*Qs0*abs(Qs0) + r01*Q01*abs(Q01) + r16*Q16*abs(Q16) + r69*Q69*abs(Q69) + r9EA*Q9EA*abs(Q9EA) - 10) / (2 * (rs0*abs(Qs0) + r01*abs(Q01) + r16*abs(Q16) + r69*abs(Q69) + r9EA*abs(Q9EA)));
		
		if (2 * (rs0*abs(Qs0) + r03*abs(Q03) + r3e*abs(Q3e) + ref*abs(Qef) + rfh*abs(Qfh) + rhEB*abs(QhEB)) > 0)
			R[10] = -(rs0*Qs0*abs(Qs0) + r03*Q03*abs(Q03) + r3e*Q3e*abs(Q3e) + ref*Qef*abs(Qef) + rfh*Qfh*abs(Qfh) + rhEB*QhEB*abs(QhEB) -10) / (2 * (rs0*abs(Qs0) + r03*abs(Q03) + r3e*abs(Q3e) + ref*abs(Qef) + rfh*abs(Qfh) + rhEB*abs(QhEB)));

		/// add the adjust flow to each Q pipes
		Q01 = Q01 + R[0] + R[9];
		Q12 = Q12 + R[0] - R[1];
		Q24 = Q24 + R[0] - R[8];
		Q34 = Q34 - R[0] + R[2];
		Q03 = Q03 - R[0] + R[10];
		Q16 = Q16 + R[1] + R[9];
		Q56 = Q56 - R[1] + R[3];
		Q25 = Q25 - R[1] + R[8];
		Q47 = Q47 - R[8] + R[2];
		Q78 = Q78 + R[2] - R[3];
		Q8a = Q8a + R[2] - R[4];

		Q57 = Q57 - R[3] + R[8];
		Q89 = Q89 - R[3] + R[4];
		Q69 = Q69 + R[3] + R[9];
		Qad = Qad - R[4] + R[5];
		Q9g = Q9g - R[4];
		Qdg = Qdg - R[4] + R[7];
		Qdc = Qdc + R[5] - R[7];
		Qgh = Qgh + R[7];
		Qfh = Qfh - R[7] + R[10];
		QhEB = QhEB + R[10];
		Q9EA = Q9EA + R[9];
		Qab = Qab + R[2] - R[5];
		Qeb = Qeb - R[2] + R[6];
		Q3e = Q3e - R[2] + R[10];
		Qef = Qef - R[6] + R[10];
		Qbc = Qbc - R[5] + R[6];
		Qcf = Qcf + R[6] - R[7];
		Qs0 = Qs0 + R[10] + R[9];

	}
	clock_t t2 = clock();

	printf("time consume : %f", t2 - t1);
	//	printf("abs = %g\n",myabs(-2));
	printf("final Qs0 = %g\n", Qs0);
	printf("final Q9-EA = %g\n", Q9EA);
	printf("final Qh-EB = %g\n", QhEB);
	
	system("pause");
	return 0;
}

