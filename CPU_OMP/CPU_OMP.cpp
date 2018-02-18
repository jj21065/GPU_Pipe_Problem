// CPU_OMP.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include<windows.h>

#define N 29
#define Pr 10.0
float rs0 = 1;

float r01 = 1;
float r16 = 5;
float r03 = 5;
float r12 = 1;
float r34 = 3;

float r24 = 1;
float r25 = 1;
float r56 = 3;
float r47 = 1;
float r3e = 1;

float r57 = 1;
float r78 = 2;
float r69 = 3;
float r89 = 3;

float r9EA = 1;

float r8a = 1;
float rab = 3;
float rad = 2;
float r9g = 1;

float reb = 2;
float rdg = 5;
float rbc = 2;
float rdc = 1;

float ref = 5;
float rcf = 2;
float rgh = 1;
float rfh = 2;
float rhEB = 1;
void Hardy_Cross_method();

void updateQ();
float R[11] = { 0 };
float* Q = new float[N];
float* r = new float[N];
int _tmain(int argc, _TCHAR* argv[])
{
	
	// pipe parameter 
	float init_Q = sqrt(Pr);
	float Q01 = init_Q / 2.0;
	float Q24 = init_Q / 8.0;
	float Q57 = init_Q / 8.0;
	float Q16 = init_Q / 4.0;
	float Q25 = init_Q / 8.0;
	float Q78 = init_Q / 4.0;
	float Q12 = init_Q / 4.0;
	float Q56 = 0;
	float Q89 = 0;
	float Q03 = init_Q / 2.0;
	float Q3e = init_Q / 2.0;
	float Q69 = init_Q / 4.0;
	float Q34 = 0;
	float Q47 = init_Q / 8.0;
	float Q8a = init_Q / 4.0;
	float Qab = init_Q / 8.0;
	float Qad = init_Q / 8.0;
	float Qeb = 0;
	float Qbc = init_Q / 8.0;
	float Qdc = 0;
	float Qdg = init_Q / 8.0;
	float Q9g = 0;
	float Qgh = init_Q / 8.0;
	float Qef = init_Q / 2.0;
	float Qcf = init_Q / 8.0;
	float Qfh = init_Q * 5.0 / 8.0;
	float QS0 = init_Q;
	float QEA = init_Q / 4.0;
	float QEB = init_Q * 3.0 / 4.0;


	Q[0] = Q01;
	Q[1] = Q24;
	Q[2] = Q57;
	Q[3] = Q16;
	Q[4] = Q25;
	Q[5] = Q78;
	Q[6] = Q12;
	Q[7] = Q56;
	Q[8] = Q89;
	Q[9] = Q03;
	Q[10] = Q3e;
	Q[11] = Q69;
	Q[12] = Q34;
	Q[13] = Q47;
	Q[14] = Q8a;
	Q[15] = Qab;
	Q[16] = Qad;
	Q[17] = Qeb;
	Q[18] = Qbc;
	Q[19] = Qdc;
	Q[20] = Qdg;
	Q[21] = Q9g;
	Q[22] = Qgh;
	Q[23] = Qef;
	Q[24] = Qcf;
	Q[25] = Qfh;
	Q[26] = QS0;
	Q[27] = QEA;
	Q[28] = QEB;

	r[0] = r01;
	r[1] = r24;
	r[2] = r57;
	r[3] = r16;
	r[4] = r25;
	r[5] = r78;
	r[6] = r12;
	r[7] = r56;
	r[8] = r89;
	r[9] = r03;
	r[10] = r3e;
	r[11] = r69;
	r[12] = r34;
	r[13] = r47;
	r[14] = r8a;
	r[15] = rab;
	r[16] = rad;
	r[17] = reb;
	r[18] = rbc;
	r[19] = rdc;
	r[20] = rdg;
	r[21] = r9g;
	r[22] = rgh;
	r[23] = ref;
	r[24] = rcf;
	r[25] = rfh;
	r[26] = rs0;
	r[27] = r9EA;
	r[28] = rhEB;

	// initial conditions 

	int iter_no = 250;
	int i;
	float Error = 1e5;

	//SYSTEMTIME t1, t2;
	//GetLocalTime(&t1);
	omp_set_num_threads(4);
	



	for (i = 0; i < iter_no; i++)
	{
#pragma omp parallel shared(R,Q,r)
	{
		Hardy_Cross_method();
	}
	updateQ();


	}



	//GetLocalTime(&t2);
	//float time = t2.wSecond - t1.wSecond + (t2.wMilliseconds - t1.wMilliseconds) / 1000.0;
	//printf("iter : %d, time : %g\n", i, time);
	for (int i = 0; i < 29; i++)
		printf("final Q[%d] = %g\n", i, Q[i]);

	system("pause");

	delete[] Q;

	return 0;
}



void updateQ()
{
	
		Q[0] = Q[0] + R[0] + R[9];
		Q[1] = Q[1] + R[0] - R[2];
		Q[2] = Q[2] + R[2] - R[4];
		Q[3] = Q[3] + R[1] + R[9];
		Q[4] = Q[4] - R[1] + R[2];
		Q[5] = Q[5] + R[3] - R[4];
		Q[6] = Q[6] + R[0] - R[1];
		Q[7] = Q[7] - R[1] + R[4];
	
		Q[8] = Q[8] - R[4] + R[5];
		Q[9] = Q[9] - R[0] + R[10];
		Q[10] = Q[10] - R[3] + R[10];
		Q[11] = Q[11] + R[4] + R[9];
		Q[12] = Q[12] - R[0] + R[3];
		Q[13] = Q[13] - R[2] + R[3];
		Q[14] = Q[14] + R[3] - R[5];
		Q[15] = Q[15] + R[3] - R[6];
		Q[16] = Q[16] + R[6] - R[5];

	
		Q[17] = Q[17] - R[3] + R[7];
		Q[18] = Q[18] - R[6] + R[7];
		Q[19] = Q[19] + R[6] - R[8];
		Q[20] = Q[20] - R[5] + R[8];
		Q[21] = Q[21] - R[5];
		Q[22] = Q[22] + R[8];
		Q[23] = Q[23] - R[7] + R[10];
		Q[24] = Q[24] + R[7] - R[8];

		Q[25] = Q[25] - R[8] + R[10];
		Q[26] = Q[26] + R[9] + R[10];
		Q[27] = Q[27] + R[9];
		Q[28] = Q[28] + R[10];
	
}
void Hardy_Cross_method()
{
	int tid = omp_get_thread_num();

	if (tid == 0)
	{
		R[0] = -(r[0] * Q[0] * abs(Q[0]) + r[6] * Q[6] * abs(Q[6]) + r[1] * Q[1] * abs(Q[1]) - r[12] * Q[12] * abs(Q[12]) - r[9] * Q[9] * abs(Q[9]))
			/ (2 * r[0] * fabs(Q[0]) + 2 * r[6] * fabs(Q[6]) + 2 * r[1] * fabs(Q[1]) + 2 * r[12] * fabs(Q[12]) + 2 * r[9] * fabs(Q[9]));

		R[1] = -(r[3] * Q[3] * abs(Q[3]) - r[7] * Q[7] * abs(Q[7]) - r[4] * Q[4] * abs(Q[4]) - r[6] * Q[6] * abs(Q[6]))
			/ (2 * r[3] * fabs(Q[3]) + 2 * r[7] * fabs(Q[7]) + 2 * r[4] * fabs(Q[4]) + 2 * r[6] * fabs(Q[6]));

		R[2] = -(r[4] * Q[4] * abs(Q[4]) + r[2] * Q[2] * abs(Q[2]) - r[13] * Q[13] * abs(Q[13]) - r[1] * Q[1] * abs(Q[1]))
			/ (2 * r[4] * fabs(Q[4]) + 2 * r[2] * fabs(Q[2]) + 2 * r[13] * fabs(Q[13]) + 2 * r[1] * fabs(Q[1]));
	}
	else if (tid == 1)
	{
		R[3] = -(r[12] * Q[12] * abs(Q[12]) + r[13] * Q[13] * abs(Q[13]) + r[5] * Q[5] * abs(Q[5]) + r[14] * Q[14] * abs(Q[14]) + r[15] * Q[15] * abs(Q[15]) - r[17] * Q[17] * abs(Q[17]) - r[10] * Q[10] * abs(Q[10]))
			/ (2 * r[12] * fabs(Q[12]) + 2 * r[13] * fabs(Q[13]) + 2 * r[5] * fabs(Q[5]) + 2 * r[14] * fabs(Q[14]) + 2 * r[15] * fabs(Q[15]) + 2 * r[17] * fabs(Q[17]) + 2 * r[10] * fabs(Q[10]));

		R[4] = -(r[7] * Q[7] * abs(Q[7]) + r[11] * Q[11] * abs(Q[11]) - r[8] * Q[8] * abs(Q[8]) - r[5] * Q[5] * abs(Q[5]) - r[2] * Q[2] * abs(Q[2]))
			/ (2 * r[7] * fabs(Q[7]) + 2 * r[11] * fabs(Q[11]) + 2 * r[8] * fabs(Q[8]) + 2 * r[5] * fabs(Q[5]) + 2 * r[2] * fabs(Q[2]));

		R[5] = -(r[8] * Q[8] * abs(Q[8]) - r[21] * Q[21] * abs(Q[21]) - r[20] * Q[20] * abs(Q[20]) - r[16] * Q[16] * abs(Q[16]) - r[14] * Q[14] * abs(Q[14]))
			/ (2 * r[8] * fabs(Q[8]) + 2 * r[21] * fabs(Q[21]) + 2 * r[20] * fabs(Q[20]) + 2 * r[16] * fabs(Q[16]) + 2 * r[14] * fabs(Q[14]));
	}

	else if (tid == 2)
	{
		R[6] = -(r[16] * Q[16] * abs(Q[16]) + r[19] * Q[19] * abs(Q[19]) - r[18] * Q[18] * abs(Q[18]) - r[15] * Q[15] * abs(Q[15]))
			/ (2 * r[16] * fabs(Q[16]) + 2 * r[19] * fabs(Q[19]) + 2 * r[18] * fabs(Q[18]) + 2 * r[15] * fabs(Q[15]));

		R[7] = -(r[17] * Q[17] * abs(Q[17]) + r[18] * Q[18] * abs(Q[18]) + r[24] * Q[24] * abs(Q[24]) - r[23] * Q[23] * abs(Q[23]))
			/ (2 * r[17] * fabs(Q[17]) + 2 * r[18] * fabs(Q[18]) + 2 * r[24] * fabs(Q[24]) + 2 * r[23] * fabs(Q[23]));

		R[8] = -(r[20] * Q[20] * abs(Q[20]) + r[22] * Q[22] * abs(Q[22]) - r[25] * Q[25] * abs(Q[25]) - r[24] * Q[24] * abs(Q[24]) - r[19] * Q[19] * abs(Q[19]))
			/ (2 * r[20] * fabs(Q[20]) + 2 * r[22] * fabs(Q[22]) + 2 * r[25] * fabs(Q[25]) + 2 * r[24] * fabs(Q[24]) + 2 * r[19] * fabs(Q[19]));

	}
	else if (tid == 3)
	{
		R[9] = -(-Pr + r[26] * Q[26] * abs(Q[26]) + r[0] * Q[0] * abs(Q[0]) + r[3] * Q[3] * abs(Q[3]) + r[11] * Q[11] * abs(Q[11]) + r[27] * Q[27] * abs(Q[27]))
			/ (2 * r[26] * fabs(Q[26]) + 2 * r[0] * fabs(Q[0]) + 2 * r[3] * fabs(Q[3]) + 2 * r[11] * fabs(Q[11]) + 2 * r[27] * fabs(Q[27]));

		R[10] = -(-Pr + r[26] * Q[26] * abs(Q[26]) + r[9] * Q[9] * abs(Q[9]) + r[10] * Q[10] * abs(Q[10]) + r[23] * Q[23] * abs(Q[23]) + r[25] * Q[25] * abs(Q[25]) + r[28] * Q[28] * abs(Q[28]))
			/ (2 * r[26] * fabs(Q[26]) + 2 * r[9] * fabs(Q[9]) + 2 * r[10] * fabs(Q[10]) + 2 * r[23] * fabs(Q[23]) + 2 * r[25] * fabs(Q[25]) + 2 * r[28] * fabs(Q[28]));

	}

//#pragma omp barrier

	
//#pragma omp barrier

	

}
