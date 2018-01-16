
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include<windows.h>

#define N 29
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
float r57 = 1;
float r78 = 2;
float r69 = 3;
float r89 = 3;
float r9EA = 1;
float r3e = 1;
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


void MatrixVectorProduct(float *a, float*p, float *ap, int n);

float Dot(float *v1, float *v2, int n);

void Sum_Store(float *store, float *x, float scalar, float*v, int n);

void conjgrad(float*A, float*x, float*b);

void Newton_Raphson_method(float *x, float*r);

void Compute_R(float*x, float*r, float *R);

void Compute_J(float*x, float*r, int n);

void Compute_invJ_mul_R(float*x, float*r, float*R, float*invJ_R);

void inverse(float*A, float* invA, int n);

float Hardy_Cross_method(float* Q);

int main()
{
	float* Q = new float[N];
	float* r = new float[N];
	// pipe parameter 
	Q[0] = sqrt(10.0 / rs0);
	Q[1] = Q[0];
	Q[2] = Q[0];
	Q[3] = 0;
	Q[4] = 0;
	Q[5] = 0;
	Q[6] = 0;
	Q[7] = 0;
	Q[8] = 0;
	Q[9] = 0;
	Q[10] = 0;

	Q[11] = 0;
	Q[12] = 0;
	Q[13] = 0;
	Q[14] = Q[0];

	Q[15] = 0;
	Q[16] = 0;
	Q[17] = 0;
	Q[18] = 0;

	Q[19] = 0;
	Q[20] = 0;
	Q[21] = 0;
	Q[22] = 0;

	Q[23] = 0;
	Q[24] = 0;
	Q[25] = 0;
	Q[26] = 0;
	
	Q[27] = Q[0];
	Q[28] = 0;

	/*Q[0] = 1;
	Q[1] = 1;
	Q[2] = 1;
	Q[3] = 1;
	Q[4] = 1;
	Q[5] = 1;
	Q[6] = 1;
	Q[7] = 1;
	Q[8] = 1;
	Q[9] = 1;
	Q[10] = 1;

	Q[11] = 1;
	Q[12] = 1;
	Q[13] = 1;
	Q[14] = 1;

	Q[15] = 1;
	Q[16] =1;
	Q[17] = 1;
	Q[18] =1;

	Q[19] = 1;
	Q[20] = 1;
	Q[21] = 1;
	Q[22] = 1;

	Q[23] = 1;
	Q[24] = 1;
	Q[25] = 1;
	Q[26] = 1;

	Q[27] = 1;
	Q[28] = 1;*/

	r[0] = 1;
	r[1] = 1;
	r[2] = 5;
	r[3] = 5;
	r[4] = 1;

	r[5] = 3;
	r[6] = 1;
	r[7] = 1;
	r[8] = 3;
	r[9] = 1;
	
	r[10] = 1;
	r[11] = 2;
	r[12] = 3;
	r[13] = 3;
	r[14] = 1;

	r[15] = 1;
	r[16] = 1;
	r[17] = 3;
	r[18] = 2;
	r[19] = 1;

	r[20] = 2;
	r[21] = 5;
	r[22] = 2;
	r[23] = 1;

	r[24] = 5;
	r[25] = 2;
	r[26] = 1;
	r[27] = 2;
	r[28] = 1;

	// initial conditions 

	int iter_no = 1000;
	int i;
	float Error = 1e5;

	SYSTEMTIME t1, t2;
	GetLocalTime(&t1);

	for (i = 0; i < iter_no; i++)
	{
		Error = Hardy_Cross_method(Q);
		/*if (Error < 1e-8)
			break;*/
		//Newton_Raphson_method(Q, r);
	}
	GetLocalTime(&t2);
	float time = t2.wSecond - t1.wSecond + (t2.wMilliseconds - t1.wMilliseconds) / 1000.0;
	printf("iter : %d, time : %g\n", i, time);
	for (int i = 0; i < 29; i++)
		printf("final Q[%d] = %g\n", i, Q[i]);

	/*
		printf("final Q[0] = %g\n", Q[0]);
		printf("final Q9-EA = %g\n", Q[27]);
		printf("final Qh-EB = %g\n", Q[28]);*/
	system("pause");

	delete[] Q;
	delete[] r;

	return 0;
}

float Hardy_Cross_method(float* Q)
{
	float R[11] = { 0 };

	/// calculate adjust flow Q
	if (2 * (r01*abs(Q[1]) + r12*abs(Q[3]) + r24*abs(Q[6]) + r34*abs(Q[5]) + r03*abs(Q[4])) > 0)
		R[0] = -(r01*Q[1] * abs(Q[1]) + r12*Q[3] * abs(Q[3]) + r24*Q[6] * abs(Q[6]) - r34*Q[5] * abs(Q[5]) - r03*Q[4] * abs(Q[4])) / (2 * (r01*abs(Q[1]) + r12*abs(Q[3]) + r24*abs(Q[6]) + r34*abs(Q[5]) + r03*abs(Q[4])));

	if (2 * (r16*abs(Q[2]) + r56*abs(Q[8]) + r25*abs(Q[7]) + r12*abs(Q[3])) > 0)
		R[1] = -(r16*Q[2] * abs(Q[2]) - r56*Q[8] * abs(Q[8]) - r25*Q[7] * abs(Q[7]) - r12*Q[3] * abs(Q[3])) / (2 * (r16*abs(Q[2]) + r56*abs(Q[8]) + r25*abs(Q[7]) + r12*abs(Q[3])));

	if (2 * (r34*abs(Q[5]) + r47*abs(Q[10]) + r78*abs(Q[12]) + r8a*abs(Q[15]) + rab*abs(Q[19]) + reb*abs(Q[16]) + r3e*abs(Q[9])) > 0)
		R[2] = -(r34*Q[5] * abs(Q[5]) + r47*Q[10] * abs(Q[10]) + r78*Q[12] * abs(Q[12]) + r8a*Q[15] * abs(Q[15]) + rab*Q[19] * abs(Q[19]) - reb*Q[16] * abs(Q[16]) - r3e*Q[9] * abs(Q[9])) / (2 * (r34*abs(Q[5]) + r47*abs(Q[10]) + r78*abs(Q[12]) + r8a*abs(Q[15]) + rab*abs(Q[19]) + reb*abs(Q[16]) + r3e*abs(Q[9])));

	if (2 * (r56*abs(Q[8]) + r69*abs(Q[14]) + r89*abs(Q[13]) + r78*abs(Q[12]) + r57*abs(Q[11])) > 0)
		R[3] = -(r56*Q[8] * abs(Q[8]) + r69*Q[14] * abs(Q[14]) - r89*Q[13] * abs(Q[13]) - r78*Q[12] * abs(Q[12]) - r57*Q[11] * abs(Q[11])) / (2 * (r56*abs(Q[8]) + r69*abs(Q[14]) + r89*abs(Q[13]) + r78*abs(Q[12]) + r57*abs(Q[11])));

	if (2 * (r89*abs(Q[13]) + r9g*abs(Q[21]) + rdg*abs(Q[17]) + rad*abs(Q[23]) + r8a*abs(Q[15])) > 0)
		R[4] = -(r89*Q[13] * abs(Q[13]) - r9g*Q[21] * abs(Q[21]) - rdg*Q[17] * abs(Q[17]) - rad*Q[23] * abs(Q[23]) - r8a*Q[15] * abs(Q[15])) / (2 * (r89*abs(Q[13]) + r9g*abs(Q[21]) + rdg*abs(Q[17]) + rad*abs(Q[23]) + r8a*abs(Q[15])));

	if (2 * (rad*abs(Q[23]) + rdc*abs(Q[24]) + rbc*abs(Q[20]) + rab*abs(Q[19])) > 0)
		R[5] = -(rad*Q[23] * abs(Q[23]) + rdc*Q[24] * abs(Q[24]) - rbc*Q[20] * abs(Q[20]) - rab*Q[19] * abs(Q[19])) / (2 * (rad*abs(Q[23]) + rdc*abs(Q[24]) + rbc*abs(Q[20]) + rab*abs(Q[19])));

	if (2 * (reb*abs(Q[16]) + rbc*abs(Q[20]) + rcf*abs(Q[22]) + ref*abs(Q[18])) > 0)
		R[6] = -(reb*Q[16] * abs(Q[16]) + rbc*Q[20] * abs(Q[20]) + rcf*Q[22] * abs(Q[22]) - ref*Q[18] * abs(Q[18])) / (2 * (reb*abs(Q[16]) + rbc*abs(Q[20]) + rcf*abs(Q[22]) + ref*abs(Q[18])));

	if (2 * (rdg*abs(Q[17]) + rgh*abs(Q[25]) + rfh*abs(Q[26]) + rcf*abs(Q[22]) + rdc*abs(Q[24])) > 0)
		R[7] = -(rdg*Q[17] * abs(Q[17]) + rgh*Q[25] * abs(Q[25]) - rfh*Q[26] * abs(Q[26]) - rcf*Q[22] * abs(Q[22]) - rdc*Q[24] * abs(Q[24])) / (2 * (rdg*abs(Q[17]) + rgh*abs(Q[25]) + rfh*abs(Q[26]) + rcf*abs(Q[22]) + rdc*abs(Q[24])));

	if (2 * (r25*abs(Q[7]) + r57*abs(Q[11]) + r47*abs(Q[10]) + r24*abs(Q[6])) > 0)
		R[8] = -(r25*Q[7] * abs(Q[7]) + r57*Q[11] * abs(Q[11]) - r47*Q[10] * abs(Q[10]) - r24*Q[6] * abs(Q[6])) / (2 * (r25*abs(Q[7]) + r57*abs(Q[11]) + r47*abs(Q[10]) + r24*abs(Q[6])));

	if (2 * (rs0*abs(Q[0]) + r01*abs(Q[1]) + r16*abs(Q[2]) + r69*abs(Q[14]) + r9EA*abs(Q[27])) > 0)
		R[9] = -(rs0*Q[0] * abs(Q[0]) + r01*Q[1] * abs(Q[1]) + r16*Q[2] * abs(Q[2]) + r69*Q[14] * abs(Q[14]) + r9EA*Q[27] * abs(Q[27]) - 10) / (2 * (rs0*abs(Q[0]) + r01*abs(Q[1]) + r16*abs(Q[2]) + r69*abs(Q[14]) + r9EA*abs(Q[27])));

	if (2 * (rs0*abs(Q[0]) + r03*abs(Q[4]) + r3e*abs(Q[9]) + ref*abs(Q[18]) + rfh*abs(Q[26]) + rhEB*abs(Q[28])) > 0)
		R[10] = -(rs0*Q[0] * abs(Q[0]) + r03*Q[4] * abs(Q[4]) + r3e*Q[9] * abs(Q[9]) + ref*Q[18] * abs(Q[18]) + rfh*Q[26] * abs(Q[26]) + rhEB*Q[28] * abs(Q[28]) - 10) / (2 * (rs0*abs(Q[0]) + r03*abs(Q[4]) + r3e*abs(Q[9]) + ref*abs(Q[18]) + rfh*abs(Q[26]) + rhEB*abs(Q[28])));

	/// add the adjust flow to each Q pipes
	Q[1] = Q[1] + R[0] + R[9];
	Q[3] = Q[3] + R[0] - R[1];
	Q[6] = Q[6] + R[0] - R[8];
	Q[5] = Q[5] - R[0] + R[2];
	Q[4] = Q[4] - R[0] + R[10];
	Q[2] = Q[2] + R[1] + R[9];
	Q[8] = Q[8] - R[1] + R[3];
	Q[7] = Q[7] - R[1] + R[8];
	Q[10] = Q[10] - R[8] + R[2];
	Q[12] = Q[12] + R[2] - R[3];
	Q[15] = Q[15] + R[2] - R[4];

	Q[11] = Q[11] - R[3] + R[8];
	Q[13] = Q[13] - R[3] + R[4];
	Q[14] = Q[14] + R[3] + R[9];
	Q[23] = Q[23] - R[4] + R[5];
	Q[21] = Q[21] - R[4];
	Q[17] = Q[17] - R[4] + R[7];
	Q[24] = Q[24] + R[5] - R[7];
	Q[25] = Q[25] + R[7];
	Q[26] = Q[26] - R[7] + R[10];
	Q[28] = Q[28] + R[10];
	Q[27] = Q[27] + R[9];
	Q[19] = Q[19] + R[2] - R[5];
	Q[16] = Q[16] - R[2] + R[6];
	Q[9] = Q[9] - R[2] + R[10];
	Q[18] = Q[18] - R[6] + R[10];
	Q[20] = Q[20] - R[5] + R[6];
	Q[22] = Q[22] + R[6] - R[7];
	Q[0] = Q[0] + R[10] + R[9];

	float tmpErr = 0;
	for (int i = 0; i < 11; i++)
	{
		tmpErr = tmpErr + R[i] * R[i];
	}
	return tmpErr;

}

void Newton_Raphson_method(float *x, float*r)
{

	/// Computet R 
	float invJR[N] = { 0 };
	float R[N] = { 0 };
	Compute_R(x, r, R);
	/// Compute J 
	Compute_invJ_mul_R(x, r, R, invJR);

	for (int i = 0; i < N; i++)
	{
		x[i] = x[i] - invJR[i];
	}
	// x = x - inv(J')*R'
}

void Compute_R(float*x, float*r, float *R)
{

	R[0] = (r01*x[1] * abs(x[1]) + r12*x[4] * abs(x[4]) + r24*x[6] * abs(x[6]) - r34*x[5] * abs(x[5]) - r03*x[3] * abs(x[3]));

	R[1] = (r16*x[2] * abs(x[2]) - r56*x[8] * abs(x[8]) - r25*x[7] * abs(x[7]) - r12*x[4] * abs(x[4]));

	R[2] = (r25*x[7] * abs(x[7]) + r57*x[10] * abs(x[10]) - r47*x[9] * abs(x[9]) - r24*x[6] * abs(x[6]));

	R[3] = (r34*x[5] * abs(x[5]) + r47*x[9] * abs(x[9]) + r78*x[11] * abs(x[11]) + r8a*x[16] * abs(x[16]) + rab*x[17] * abs(x[17]) - reb*x[20] * abs(x[20]) - r3e*x[15] * abs(x[15]));

	R[4] = (r56*x[8] * abs(x[8]) + r69*x[12] * abs(x[12]) - r89*x[13] * abs(x[13]) - r78*x[11] * abs(x[11]) - r57*x[10] * abs(x[10]));

	R[5] = (r89*x[13] * abs(x[13]) - r9g*x[19] * abs(x[19]) - rdg*x[21] * abs(x[21]) - rad*x[18] * abs(x[18]) - r8a*x[16] * abs(x[16]));

	R[6] = (rad*x[18] * abs(x[18]) + rdc*x[23] * abs(x[23]) - rbc*x[22] * abs(x[22]) - rab*x[17] * abs(x[17]));

	R[7] = (reb*x[20] * abs(x[20]) + rbc*x[22] * abs(x[22]) + rcf*x[25] * abs(x[25]) - ref*x[24] * abs(x[24]));

	R[8] = (rdg*x[21] * abs(x[21]) + rgh*x[26] * abs(x[26]) - rfh*x[27] * abs(x[27]) - rcf*x[25] * abs(x[25]) - rdc*x[23] * abs(x[23]));

	R[9] = x[25] + x[24] - x[27];

	R[10] = (rs0*x[0] * abs(x[0]) + r01*x[1] * abs(x[1]) + r16*x[2] * abs(x[2]) + r69*x[12] * abs(x[12]) + r9EA*x[14] * abs(x[14]) - 10);

	R[11] = (rs0*x[0] * abs(x[0]) + r03*x[3] * abs(x[3]) + r3e*x[15] * abs(x[15]) + ref*x[24] * abs(x[24]) + rfh*x[27] * abs(x[27]) + rhEB*x[28] * abs(x[28]) - 10);

	R[12] = x[0] - x[1] - x[3];

	R[13] = x[1] - x[4] - x[2];

	R[14] = x[4] - x[6] - x[7];

	R[15] = x[3] - x[5] - x[15];

	R[16] = x[5] + x[6] - x[9];

	R[17] = x[7] - x[8] - x[10];

	R[18] = x[2] + x[8] - x[12];

	R[19] = x[10] + x[9] - x[11];

	R[20] = x[11] - x[13] - x[16];

	R[21] = x[13] + x[12] + x[19] - x[14];

	R[22] = x[16] - x[17] - x[18];

	R[23] = x[15] - x[20] - x[24];

	R[24] = x[20] + x[17] - x[22];

	R[25] = x[18] - x[21] - x[23];

	R[26] = x[21] + x[19] - x[26];

	R[27] = x[26] + x[27] - x[28];

	R[28] = x[22] + x[23] - x[25];


}

void Compute_invJ_mul_R(float*x, float*r, float*R, float*invJ_R)
{
	/// Conjugate gradient method 
	float A[N*N] = { 0 };

	A[0] = 0; A[1] = 2 * r[1] * x[1];  A[4] = 2 * r[4] * x[4]; A[6] = 2 * r[6] * x[6]; A[5] = -2 * r[5] * x[5]; A[3] = -2 * r[3] * x[3];

	A[29] = 0; A[29 + 2] = 2 * r[2] * x[2]; A[29 + 8] = -2 * r[8] * x[8]; A[29 + 7] = -2 * r[7] * x[7]; A[29 + 4] = -2 * r[4] * x[4];

	A[58] = 0; A[58 + 7] = 2 * r[7] * x[7]; A[58 + 10] = 2 * r[10] * x[10]; A[58 + 9] = -2 * r[9] * x[9]; A[58 + 6] = -2 * r[6] * x[6];

	A[87] = 0; A[87 + 5] = 2 * r[5] * x[5]; A[9] = 2 * r[9] * x[9]; A[87 + 11] = 2 * r[11] * x[11]; A[87 + 16] = 2 * r[16] * x[16]; A[87 + 17] = 2 * r[17] * x[17]; A[87 + 20] = -2 * r[20] * x[20]; A[87 + 15] = -2 * r[15] * x[15];

	A[116] = 0; A[116 + 8] = 2 * r[8] * x[8]; A[116 + 12] = 2 * r[12] * x[12]; A[116 + 13] = -2 * r[13] * x[13]; A[116 + 11] = -2 * r[11] * x[11]; A[116 + 10] = -2 * r[10] * x[10];
	
	///
	
	A[145] = 0; A[145 + 13] = 2 * r[13] * x[13]; A[145 + 19] = -2 * r[19] * x[19]; A[145 + 21] = -2 * r[21] * x[21]; A[145 + 18] = -2 * r[18] * x[18]; A[145 + 16] = -2 * r[16] * x[16];

	A[174] = 0; A[174 + 18] = 2 * r[18] * x[18]; A[174 + 23] = 2 * r[23] * x[23]; A[174 + 22] = -2 * r[22] * x[22]; A[174 + 17] = -2 * r[17] * x[17];

	A[203] = 0; A[203 + 20] = 2 * r[20] * x[20]; A[203 + 22] = 2 * r[22] * x[22]; A[203 + 25] = 2 * r[25] * x[25]; A[203 + 24] = -2 * r[24] * x[24];

	A[232] = 0; A[232 + 21] = 2 * r[21] * x[21]; A[232 + 26] = 2 * r[26] * x[26]; A[232 + 27] = -2 * r[27] * x[27]; A[232 + 25] = -2 * r[25] * x[25]; A[232 + 23] = -2 * r[23] * x[23];

	A[261] = 0; A[261 + 25] = 1; A[261 + 24] = 1; A[261 + 27] = -1;

	///

	A[290] = 2 * r[0] * x[0]; A[290 + 1] = 2 * r[1] * x[1]; A[290 + 2] = 2 * r[2] * x[2]; A[290 + 12] = 2 * r[12] * x[12]; A[290 + 14] = 2 * r[14] * x[14];

	A[319] = 2 * r[0] * x[0]; A[319 + 3] = 2 * r[3] * x[3]; A[319 + 15] = 2 * r[15] * x[15]; A[319 + 24] = 2 * r[24] * x[24]; A[319 + 27] = 2 * r[27] * x[27]; A[319 + 28] = 2 * r[28] * x[28];

	A[348] = 1; A[348 + 1] = -1; A[348 + 3] = -1;

	A[377] = 0; A[377 + 1] = 1; A[377 + 4] = -1; A[377 + 2] = 1;

	A[406] = 0; A[406+4] = 1; A[406+6] = -1; A[406+7] = -1;

	A[435] = 0; A[435+3] = 1; A[435+5] = -1; A[435+15] = -1;

	A[464] = 0; A[464+5] = 1; A[464+6] = 1; A[464+9] = -1;

	A[493] = 0; A[493+7] = 1; A[493+8] = 1; A[493+10] = -1;

	A[522] = 0; A[522+2] = 1; A[522+8] = 1; A[522+12] = -1;

	A[551] = 0; A[551+10] = 1; A[551+9] = 1; A[551+11] = -1;

	A[580] = 0; A[580+11] = 1; A[580+13] = -1; A[580+16] = -1;

	A[609] = 0; A[609+13] = 1; A[609+12] = 1; A[609+19] = 1; A[609+14] = -1;

	A[638] = 0; A[638+16] = 1; A[638+17] = -1; A[638+18] = -1;

	A[667] = 0; A[667+15] = 1; A[667+20] = -1; A[667+24] = -1;

	A[696] = 0; A[696+20] = 1; A[696+17] = 1; A[696+22] = -1;

	A[725] = 0; A[725+18] = 1; A[725+21] = -1; A[725+23] = -1;

	A[754] = 0; A[754+21] = 1; A[754+19] = 1; A[754+26] = -1;

	A[783] = 0; A[783+26] = 1; A[783+27] = 1; A[783+28] = -1;

	A[812] = 0; A[812+22] = 1; A[812+23] = 1; A[812+25] = -1;

	//FILE*pfile;
	//pfile = fopen("output.txt", "w");

	//for (int i = 0; i < N; i++)
	//{
	//	for (int j = 0; j < N; j++)
	//	{
	//		fprintf(pfile, "%g\t", A[i*N + j]);
	//	}
	//	fprintf(pfile, "\n");
	//}
	//fclose(pfile);
	//	conjgrad(A, invJ_R, R);
	float Atemp[N*N] = { 0 };
	//inverse(A, Atemp, N);
	conjgrad(A, invJ_R, R);
	//MatrixVectorProduct(Atemp, R, invJ_R, N);

}

void MatrixVectorProduct(float *a, float*p, float *ap, int n)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < n; i++){
		ap[i] = 0.0;
		for (j = 0; j < n; j++){
			ap[i] = ap[i] + a[i*n + j] * p[j];
		}
	}
}

float Dot(float *v1, float *v2, int n)
{
	int i;
	float ans = 0.0;
	for (i = 0; i < n; i++){
		ans += v1[i] * v2[i];

	}
	return ans;
}

void Sum_Store(float *store, float *x, float scalar, float*v, int n)
{
	int i = 0;
	for (i = 0; i < n; i++){
		store[i] = x[i] + scalar * v[i];
	}
}

void conjgrad(float*A, float*x, float*b)
{
	float ap[N];
	float p[N];
	float rr[N];
	float rsold = 0;
	float rsnew = 0;
	int i = 0;
	MatrixVectorProduct(A, x, ap, N);
	for (i = 0; i < N; i++){
		rr[i] = b[i] - ap[i];
		p[i] = rr[i];
	}
	rsold = Dot(rr, rr, N);

	for (i = 0; i < N; i++){
		MatrixVectorProduct(A, p, ap, N);
		/*	for (int j = 0; j < N*N; j++){
				printf("ap[%d] = %g\n", j, ap[j]);
				}*/
		float tmpap = Dot(ap, p, N);
		float alpha = rsold / (tmpap);
		Sum_Store(x, x, alpha, p, N);
		Sum_Store(rr, rr, -alpha, ap, N);
		rsnew = Dot(rr, rr, N);
		if (sqrt(rsnew) < 1e-10)
			break;
		Sum_Store(p, rr, (rsnew / rsold), p, N);
		rsold = rsnew;

	}
	// 	for (i = 0; i < N; i++)
	// 		printf("x[%d] = %g\n", i, x[i]);


}

void inverse(float*A, float* invA, int n)
{
	int i = 0, j = 0, k = 0;
	float d;
	float *a = new float[n*n * 2];

	for (i = 0; i < n*n * 2; i++)
		a[i] = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a[i*n * 2 + j] = A[i*n + j];
		}
	}
	//for (i = 0; i < n; i++)
	//{
	//	for (j = 0; j < 2 * n; j++)
	//	{
	//		printf("%g\t", a[i*n*2 + j]);
	//	}
	//	printf("\n");
	//}

	for (i = 0; i < n; i++)
		for (j = 0; j < 2 * n; j++)
			if (j == (i + n))
				a[i*n * 2 + j] = 1;

	//for (i = 0; i < n; i++)
	//{
	//	for (j = 0; j < 2 * n; j++)
	//	{
	//		printf("%g\t", a[i*n*2 + j]);
	//	}
	//	printf("\n");
	//}

	/************** partial pivoting **************/
	for (i = n - 1; i > 0; i--)
	{
		if (a[(i - 1)*n * 2 + 0] < a[i*n * 2 + 0])
			for (j = 0; j < n * 2; j++)
			{
				d = a[i*n * 2 + j];
				a[i*n * 2 + j] = a[(i - 1)*n * 2 + j];
				a[(i - 1)*n * 2 + j] = d;
			}
	}
	//cout << "pivoted output: " << endl;
	/*for (i = 0; i < n; i++)
	{
	for (j = 0; j < 2 * n; j++)
	{
	printf("%g\t", a[i*n*2 + j]);
	}
	printf("\n");
	}*/
	/********** reducing to diagonal  matrix ***********/

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			if (j != i)
			{
				d = a[j*n * 2 + i] / a[i*n * 2 + i];
				for (k = 0; k < n * 2; k++)
					a[j*n * 2 + k] -= a[i*n * 2 + k] * d;
			}
	}
	//for (i = 0; i < n; i++)
	//{
	//	for (j = 0; j < 2 * n; j++)
	//	{
	//		printf("%g\t", a[i*n * 2 + j]);
	//	}
	//	printf("\n");
	//}
	/************** reducing to unit matrix *************/
	for (i = 0; i < n; i++)
	{
		d = a[i*n * 2 + i];
		if (d != 0)
			for (j = 0; j < n * 2; j++)
				a[i*n * 2 + j] = a[i*n * 2 + j] / d;
	}

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			invA[i*n + j] = a[i*n * 2 + n + j];

	}
	//for (i = 0; i < n; i++)
	//{
	//	for (j = 0; j < 2 * n; j++)
	//	{
	//		printf("%g\t", a[i*n*2 + j]);
	//	}
	//	printf("\n");
	//}

	//for (i = 0; i < n; i++)
	//{
	//	for (j = 0; j <  n; j++)
	//	{
	//		printf("%g\t", invA[i*n+j]);
	//	}
	//	printf("\n");
	//}
	delete[] a;
}

float norm(float *a, int n)
{

	float value = 0;
	for (int i = 0; i < n; i++)
	{
		value = value + a[i] * a[i];
	}
	return value;
}

//void BICG(float *A, float *x, float *b)
//{
//	int iter = N;
//	int flag = 0;
//
//	float norm2 = norm(b, N);
//	if (norm2 == 0.0)
//		norm2 = 1.0;
//	float ap[N] = { 0 };
//	float r[N] = { 0 };
//	MatrixVectorProduct(A, x, ap, N);
//	for (int i = 0; i < N; i++){
//		r[i] = b[i] - ap[i];
//	}
//		
//	float error = norm(r, N) / norm2;
//	if (error < 0.001) 
//		return;
//
//	float r_tld[N];
//	for (int i = 0; i < N; i++)
//		r_tld[i] = r[i];
//
//	for(int i = 0;i<iter;i++)
//		                   
//
//		z = M \ r;
//	z_tld = M' \ r_tld;
//		rho = (z'*r_tld );
//		if (rho == 0.0),
//			break
//			end
//
//			if (iter > 1), % direction vectors
//				beta = rho / rho_1;
//	p = z + beta*p;
//	p_tld = z_tld + beta*p_tld;
//			else
//				p = z;
//	p_tld = z_tld;
//	end
//
//		q = A*p;                            % compute residual pair
//		q_tld = A'*p_tld;
//		alpha = rho / (p_tld'*q );
//
//		x = x + alpha*p;                    % update approximation
//		r = r - alpha*q;
//	r_tld = r_tld - alpha*q_tld;
//
//	error = norm(r) / bnrm2;          % check convergence
//		if (error <= tol), break, end
//
//			rho_1 = rho;
//
//	end
//
//		if (error <= tol), % converged
//			flag = 0;
//	elseif(rho == 0.0), % breakdown
//		flag = -1;
//		else
//			flag = 1;                           % no convergence
//			end
//
//}