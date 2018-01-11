
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


void MatrixVectorProduct(float *a, float*p, float *ap, int n);

float Dot(float *v1, float *v2, int n);

void Sum_Store(float *store, float *x, float scalar, float*v, int n);

void conjgrad(float*A, float*x, float*b);

void Newton_Raphson_method(float *x, float*r);

void Compute_R(float*x, float*r, float *R);

void Compute_J(float*x, float*r, int n);

void Compute_invJ_mul_R(float*x, float*r, float*R, float*invJ_R);

int main()
{
	float* Q = new float[N];
	float* r = new float[N];
	// pipe parameter 
	Q[0] = sqrt(10.0 / rs0);
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
	Q[16] = 1;
	Q[17] = 1;
	Q[18] = 1;
	Q[19] = 1;
	Q[20] = 1;
	Q[21] = 1;
	Q[22] = 1;
	Q[23] = 1;
	Q[24] = 1;
	Q[25] = 1;
	Q[26] = 1;
	Q[27] = 1;
	Q[28] = 1;

	r[0] = 1;
	r[1] = 1;
	r[2] = 5;
	r[3] = 1;
	r[4] = 5;
	r[5] = 3;
	r[6] = 1;
	r[7] = 1;
	r[8] = 3;
	r[9] = 1;
	r[10] = 1;
	r[11] = 1;
	r[12] = 2;
	r[13] = 3;
	r[14] = 3;
	r[15] = 1;
	r[16] = 2;
	r[17] = 5;
	r[18] = 5;
	r[19] = 3;
	r[20] = 2;
	r[21] = 1;
	r[22] = 2;
	r[23] = 2;
	r[24] = 1;
	r[25] = 1;
	r[26] = 2;
	r[27] = 1;
	r[28] = 1;


	// initial conditions 

	int iter_no = 10;
	int i;


	printf("inital Qs0 = %g\n", Q[0]);
	for (i = 0; i < iter_no; i++)
	{
		//	Hardy_Cross_method();


		Newton_Raphson_method(Q, r);
	}

	printf("final Qs0 = %g\n", Q[0]);
	printf("final Q9-EA = %g\n", Q[27]);
	printf("final Qh-EB = %g\n", Q[28]);

	free(Q);
	free(r);
	system("pause");
	return 0;
}

void Hardy_Cross_method()
{
	float R[11] = { 0 };

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
		R[10] = -(rs0*Qs0*abs(Qs0) + r03*Q03*abs(Q03) + r3e*Q3e*abs(Q3e) + ref*Qef*abs(Qef) + rfh*Qfh*abs(Qfh) + rhEB*QhEB*abs(QhEB) - 10) / (2 * (rs0*abs(Qs0) + r03*abs(Q03) + r3e*abs(Q3e) + ref*abs(Qef) + rfh*abs(Qfh) + rhEB*abs(QhEB)));

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

	R[0] = (r01*x[1] * abs(x[1]) + r12*x[3] * abs(x[3]) + r24*x[6] * abs(x[6]) - r34*x[5] * abs(x[5]) - r03*x[4] * abs(x[4]));

	R[1] = (r16*x[2] * abs(x[2]) - r56*x[8] * abs(x[8]) - r25*x[7] * abs(x[7]) - r12*x[3] * abs(x[3]));

	R[2] = (r34*x[5] * abs(x[5]) + r47*x[10] * abs(x[10]) + r78*x[12] * abs(x[12]) + r8a*x[15] * abs(x[15]) + rab*x[19] * abs(x[19]) - reb*x[16] * abs(x[16]) - r3e*x[9] * abs(x[9]));

	R[3] = (r56*x[8] * abs(x[8]) + r69*x[14] * abs(x[14]) - r89*x[13] * abs(x[13]) - r78*x[12] * abs(x[12]) - r57*x[11] * abs(x[11]));

	R[4] = (r89*x[13] * abs(x[13]) - r9g*x[21] * abs(x[21]) - rdg*x[17] * abs(x[17]) - rad*x[23] * abs(x[23]) - r8a*x[15] * abs(x[15]));

	R[5] = (rad*x[23] * abs(x[23]) + rdc*x[24] * abs(x[24]) - rbc*x[20] * abs(x[20]) - rab*x[19] * abs(x[19]));

	R[6] = (reb*x[16] * abs(x[16]) + rbc*x[20] * abs(x[20]) + rcf*x[22] * abs(x[22]) - ref*x[18] * abs(x[18]));

	R[7] = (rdg*x[17] * abs(x[17]) + rgh*x[25] * abs(x[25]) - rfh*x[26] * abs(x[26]) - rcf*x[22] * abs(x[22]) - rdc*x[24] * abs(x[24]));

	R[8] = (r25*x[7] * abs(x[7]) + r57*x[11] * abs(x[11]) - r47*x[10] * abs(x[10]) - r24*x[6] * abs(x[6]));

	R[9] = (rs0*x[0] * abs(x[0]) + r01*x[1] * abs(x[1]) + r16*x[2] * abs(x[2]) + r69*x[14] * abs(x[14]) + r9EA*x[27] * abs(x[27]) - 10);

	R[10] = (rs0*x[0] * abs(x[0]) + r03*x[4] * abs(x[4]) + r3e*x[9] * abs(x[9]) + ref*x[18] * abs(x[18]) + rfh*x[26] * abs(x[26]) + rhEB*x[28] * abs(x[28]) - 10);

	R[11] = x[1] + x[4] - x[0];
	
	R[12] = x[2] + x[3] - x[1];
	
	R[14] = x[6] + x[7] - x[3];

	R[15] = x[5] + x[9] - x[4];

	R[16] = x[5] + x[6] - x[10];

	R[17] = x[10] + x[11] - x[12];

	R[18] = x[13] + x[15] - x[12];

	R[19] = x[8] + x[11] - x[7];

	R[20] = x[2] + x[8] - x[14];

	R[21] = x[14] + x[13] + x[21] - x[27];

	R[22] = x[19] + x[23] - x[15];

	R[23] = x[16] + x[19] - x[20];

	R[24] = x[16] + x[19] - x[9];

	R[25] = x[20] + x[24] - x[22];

	R[26] = x[22] + x[18] - x[26];

	R[27] = x[21] + x[25] - x[17];

	R[28] = x[25] + x[26] - x[28];
	

}

void Compute_invJ_mul_R(float*x, float*r,float*R,float*invJ_R)
{

	/// Conjugate gradient method 
	float A[N*N] = { 0 };

	A[0] = 0; A[1] = 2 * r[1] * x[1];  A[3] = 2 * r[3] * x[3]; A[4] = -2 * r[4] * x[4]; A[5] = -2 * r[5] * x[5]; A[6] = 2 * r[6] * x[6];
	
	A[29] = 0; A[30] = 2 * r[2] * x[2]; A[31] = -2 * r[3] * x[3]; A[36] = 2 * r[7] * x[7]; A[37] = 2 * r[8] * x[8];
	
	A[58] = 0; A[63] = 2 * r[5] * x[5]; A[67] = -2 * r[9] * x[9]; A[68] = 2 * r[10] * x[10]; A[70] = 2 * r[12] * x[12]; A[73] = 2 * r[15] * x[15];A[74] = -2*r[16]*x[16]; A[77] = 2 * r[19] * x[19];
	
	A[87] = 0; A[95] = 2 * r[8] * x[8]; A[98] = -2 * r[11] * x[11]; A[99] = -2 * r[12] * x[12]; A[100] = -2 * r[13] * x[13]; A[101] = 2 * r[14] * x[14];
	
	A[116] = 0; A[129] = 2 * r[13] * x[13]; A[131] = -2 * r[15] * x[15]; A[133] = -2 * r[17] * x[17]; A[137] = 2 * r[21] * x[21]; A[139] = 2 * r[23] * x[23];
	
	A[145] = 0; A[164] = -2 * r[19] * x[19]; A[165] = -2 * r[20] * x[20]; A[168] = 2 * r[23] * x[23]; A[169] = 2 * r[24] * x[24];

	A[174] = 0; A[190] = 2 * r[16] * x[16]; A[192] = -2 * r[18] * x[18]; A[194] = 2 * r[20] * x[20]; A[196] = 2 * r[22] * x[22];

	A[203] = 0; A[220] = 2 * r[17] * x[17]; A[225] = -2 * r[22] * x[22]; A[227] = -2 * r[24] * x[24]; A[228] = 2 * r[25] * x[25]; A[229] = -2 * r[26] * x[26];

	A[232] = 0; A[238] = -2 * r[6] * x[6]; A[239] = 2 * r[7] * x[7]; A[242] = -2 * r[10] * x[10]; A[243] = 2 * r[11] * x[11];
	
	A[261] = 2 * r[0] * x[0]; A[262] = 2 * r[1] * x[1]; A[263] = 2 * r[2] * x[2]; A[275] = 2 * r[14] * x[14]; A[288] = 2 * r[27] * x[27];

	A[290] = 2 * r[0] * x[0]; A[294] = 2 * r[4] * x[4]; A[299] = 2 * r[9] * x[9]; A[308] = 2 * r[18] * x[18]; A[316] = 2 * r[26] * x[26]; A[318] = 2 * r[28] * x[28];

	A[319] = -1; x[320] = 1; x[323] = 1;

	A[348] = 0; A[349] = -1; A[350] = 1; A[351] = 1;

	A[377] = 0; A[380] = -1; A[383] = 1; A[384] = 1;

	A[406] = 0; A[410] = -1; A[411] = 1; A[415] = 1;

	A[435] = 0; A[440] = 1; A[441] = 1; A[445] = -1; 

	A[464] = 0; A[474] = 1; A[475] = 1; A[476] = -1;

	A[493] = 0; A[503] = 1; A[504] = 1; A[505] = -1;

	A[522] = 0; A[534] = -1; A[535] = 1; A[537] = 1;

	A[551] = 0; A[558] = -1; A[559] = 1; A[562] = 1;

	A[580] = 0; A[582] = 1; A[588] = 1; A[594] = -1;

	A[609] = 0; A[622] = 1; A[623] = 1; A[630] = 1; A[636] = -1;

	A[638] = 0; A[653] = -1; A[657] = 1; A[661] = 1;

	A[667] = 0; A[683] = 1; A[686] = 1; A[687] = -1;

	A[696] = 0; A[705] = -1; A[712] = 1; A[715] = 1;

	A[725] = 0; A[745] = 1; A[749] = 1; A[747] = -1;

	A[754] = 0; A[772] = 1; A[776] = 1; A[780] = -1;

	A[783] = 0; A[800] = -1; A[804] = 1; A[808] = 1;
	
	A[812] = 0; A[837] = 1; A[838] = 1; A[840] = -1;

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
	conjgrad(A, invJ_R, R);

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