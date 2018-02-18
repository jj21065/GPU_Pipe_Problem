
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda.h>
#include <device_functions.h>
#include <cuda_runtime_api.h>
#pragma once
#ifdef __INTELLISENSE__
void __syncthreads();

#endif
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<windows.h>

#define N 11
#define NN 29
#define TPB 1
#define Pr 10.0
#define BPG (N+TPB-1)/TPB
float *Q; // 29
float *d_Q; // 29
float *dR;// 11
void Allocate_Memory(int n)
{
	size_t size = NN*sizeof(float);
	
	Q = (float*)malloc(size);
	cudaError_t error = cudaMalloc((void**)&d_Q, size);
	printf("Allocate mem dQ: %s\n", cudaGetErrorString(error));
	error = cudaMalloc((void**)&dR, N*sizeof(float));
	printf("Allocate mem dR : %s\n", cudaGetErrorString(error));
}

void Free_Memory()
{
	if (Q)
		free(Q);
	cudaError_t error = cudaFree(d_Q);
	printf("Free mem dQ : %s\n", cudaGetErrorString(error));
	error = cudaFree(dR);
	printf("Free mem dR: %s\n", cudaGetErrorString(error));
}

void CopyMemToDevice(int n)
{
	/*for (int i = 0; i < n; i++)
	{
	R[i] = data[i];
	}*/
	size_t size = NN*sizeof(float);
	cudaError_t error = cudaMemcpy(d_Q, Q, size, cudaMemcpyHostToDevice);
	printf("Memcpy Host to Device : %s\n", cudaGetErrorString(error));
}

void CopyMemToHost( int n)
{
	size_t size = NN*sizeof(float);
	cudaError_t error = cudaMemcpy(Q, d_Q, size, cudaMemcpyDeviceToHost);
	printf("Memcpy Device to Host : %s\n", cudaGetErrorString(error));
	/*for (int i = 0; i < n; i++)
	{
	data[i] = R[i];
	}*/
}

__global__ void Compute_Q(float* pQ, float*dR, int n);
__global__ void Add_Q(float* pQ, float*dR, int n);
__global__ void ErrorSum(float *error, int n);
__global__ void getError(float* dR, float *error, int n);
int main()
{
	
	
	Allocate_Memory(N);
	
	cudaMemset(dR, 0.0,N*sizeof(float));

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

	int n = 2;
	int iter_no = 100;
	int i;
	float h_error = 1e5;
	float*d_error;
	cudaError_t tmperror = cudaMalloc((void**)&d_error, BPG*sizeof(float));
	SYSTEMTIME t1,t2;
	GetLocalTime(&t1);
	
	//// main computation 

	
	CopyMemToDevice(NN);
	for (i = 0; i < 1000;i++)
	{
		Compute_Q<<<BPG,TPB>>>(d_Q,dR, N);
		getError<<<BPG,TPB>>>(dR, d_error, N);
		ErrorSum << <1, 1 >> >(d_error, BPG);
		cudaMemcpy(&h_error, &(d_error[0]), sizeof(float), cudaMemcpyDeviceToHost);
		if (h_error < 1e-12)
		{
			break;
		}
		Add_Q << <1, 1 >> >(d_Q, dR, N);
		
	}
	CopyMemToHost(NN);
	/////

	GetLocalTime(&t2);
	float time = t2.wSecond - t1.wSecond + (t2.wMilliseconds - t1.wMilliseconds) / 1000.0;
	printf("iter : %d, time : %g\n", i,time);
	//printf("time consume : %f\n", times);

	for (int i = 0; i < NN; i++)
	{
		printf("Q[%d] : %g\n", i ,Q[i]);
	}
	
	cudaFree(d_error);
	system("pause");
	Free_Memory();
	return 0;
}
__global__ void ErrorSum(float *error,int n)
{
	int i = 0;
	for (i = 0; i < n; i++)
	{
		error[0] = error[0] + error[i];
	}
}
__global__ void getError(float* dR,float *error, int n)
{
	
	__shared__ float c[TPB];
	int I = threadIdx.x;
	int i = TPB*blockIdx.x + I;
	
	if (i < n)
		c[I] = dR[i]*dR[i];
	for (int stride = blockDim.x / 2; stride > 0; stride = stride / 2)
	{
		if (I < stride){
			c[I] += c[I + stride];
		}
		__syncthreads();
	}

	if (I == 0)
		error[blockIdx.x] = c[0];
}
__global__ void Compute_Q(float* pQ, float*dR, int n)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	float r[29] = { 0 };
	// pipe parameter 
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
	//float R[11] = { 0 };

	if (i < n)
	{
			// calculate adjust flow pQ
		if (i == 0)
		{
			dR[0] = -(r[0] * pQ[0] * abs(pQ[0]) + r[6] * pQ[6] * abs(pQ[6]) + r[1] * pQ[1] * abs(pQ[1]) - r[12] * pQ[12] * abs(pQ[12]) - r[9] * pQ[9] * abs(pQ[9]))
				/ (2 * r[0] * fabs(pQ[0]) + 2 * r[6] * fabs(pQ[6]) + 2 * r[1] * fabs(pQ[1]) + 2 * r[12] * fabs(pQ[12]) + 2 * r[9] * fabs(pQ[9]));
		} 
		else if (i == 1)
		{
			dR[1] = -(r[3] * pQ[3] * abs(pQ[3]) - r[7] * pQ[7] * abs(pQ[7]) - r[4] * pQ[4] * abs(pQ[4]) - r[6] * pQ[6] * abs(pQ[6]))
				/ (2 * r[3] * fabs(pQ[3]) + 2 * r[7] * fabs(pQ[7]) + 2 * r[4] * fabs(pQ[4]) + 2 * r[6] * fabs(pQ[6]));
		}
		else if (i == 2){
			dR[2] = -(r[4] * pQ[4] * abs(pQ[4]) + r[2] * pQ[2] * abs(pQ[2]) - r[13] * pQ[13] * abs(pQ[13]) - r[1] * pQ[1] * abs(pQ[1]))
				/ (2 * r[4] * fabs(pQ[4]) + 2 * r[2] * fabs(pQ[2]) + 2 * r[13] * fabs(pQ[13]) + 2 * r[1] * fabs(pQ[1]));
		}
		else if (i == 3){
			dR[3] = -(r[12] * pQ[12] * abs(pQ[12]) + r[13] * pQ[13] * abs(pQ[13]) + r[5] * pQ[5] * abs(pQ[5]) + r[14] * pQ[14] * abs(pQ[14]) + r[15] * pQ[15] * abs(pQ[15]) - r[17] * pQ[17] * abs(pQ[17]) - r[10] * pQ[10] * abs(pQ[10]))
				/ (2 * r[12] * fabs(pQ[12]) + 2 * r[13] * fabs(pQ[13]) + 2 * r[5] * fabs(pQ[5]) + 2 * r[14] * fabs(pQ[14]) + 2 * r[15] * fabs(pQ[15]) + 2 * r[17] * fabs(pQ[17]) + 2 * r[10] * fabs(pQ[10]));
		}
		else if (i == 4){
			dR[4] = -(r[7] * pQ[7] * abs(pQ[7]) + r[11] * pQ[11] * abs(pQ[11]) - r[8] * pQ[8] * abs(pQ[8]) - r[5] * pQ[5] * abs(pQ[5]) - r[2] * pQ[2] * abs(pQ[2]))
				/ (2 * r[7] * fabs(pQ[7]) + 2 * r[11] * fabs(pQ[11]) + 2 * r[8] * fabs(pQ[8]) + 2 * r[5] * fabs(pQ[5]) + 2 * r[2] * fabs(pQ[2]));
		}
		else if (i == 5){
			dR[5] = -(r[8] * pQ[8] * abs(pQ[8]) - r[21] * pQ[21] * abs(pQ[21]) - r[20] * pQ[20] * abs(pQ[20]) - r[16] * pQ[16] * abs(pQ[16]) - r[14] * pQ[14] * abs(pQ[14]))
				/ (2 * r[8] * fabs(pQ[8]) + 2 * r[21] * fabs(pQ[21]) + 2 * r[20] * fabs(pQ[20]) + 2 * r[16] * fabs(pQ[16]) + 2 * r[14] * fabs(pQ[14]));
		}
		else if (i == 6){
			dR[6] = -(r[16] * pQ[16] * abs(pQ[16]) + r[19] * pQ[19] * abs(pQ[19]) - r[18] * pQ[18] * abs(pQ[18]) - r[15] * pQ[15] * abs(pQ[15]))
				/ (2 * r[16] * fabs(pQ[16]) + 2 * r[19] * fabs(pQ[19]) + 2 * r[18] * fabs(pQ[18]) + 2 * r[15] * fabs(pQ[15]));
		}
		else if (i == 7){
			dR[7] = -(r[17] * pQ[17] * abs(pQ[17]) + r[18] * pQ[18] * abs(pQ[18]) + r[24] * pQ[24] * abs(pQ[24]) - r[23] * pQ[23] * abs(pQ[23]))
				/ (2 * r[17] * fabs(pQ[17]) + 2 * r[18] * fabs(pQ[18]) + 2 * r[24] * fabs(pQ[24]) + 2 * r[23] * fabs(pQ[23]));
		}
		else if (i == 8)
		{
			dR[8] = -(r[20] * pQ[20] * abs(pQ[20]) + r[22] * pQ[22] * abs(pQ[22]) - r[25] * pQ[25] * abs(pQ[25]) - r[24] * pQ[24] * abs(pQ[24]) - r[19] * pQ[19] * abs(pQ[19]))
				/ (2 * r[20] * fabs(pQ[20]) + 2 * r[22] * fabs(pQ[22]) + 2 * r[25] * fabs(pQ[25]) + 2 * r[24] * fabs(pQ[24]) + 2 * r[19] * fabs(pQ[19]));
		}
		else if (i == 9){
			dR[9] = -(-Pr + r[26] * pQ[26] * abs(pQ[26]) + r[0] * pQ[0] * abs(pQ[0]) + r[3] * pQ[3] * abs(pQ[3]) + r[11] * pQ[11] * abs(pQ[11]) + r[27] * pQ[27] * abs(pQ[27]))
				/ (2 * r[26] * fabs(pQ[26]) + 2 * r[0] * fabs(pQ[0]) + 2 * r[3] * fabs(pQ[3]) + 2 * r[11] * fabs(pQ[11]) + 2 * r[27] * fabs(pQ[27]));
		}
		else if (i == 10){
			dR[10] = -(-Pr + r[26] * pQ[26] * abs(pQ[26]) + r[9] * pQ[9] * abs(pQ[9]) + r[10] * pQ[10] * abs(pQ[10]) + r[23] * pQ[23] * abs(pQ[23]) + r[25] * pQ[25] * abs(pQ[25]) + r[28] * pQ[28] * abs(pQ[28]))
				/ (2 * r[26] * fabs(pQ[26]) + 2 * r[9] * fabs(pQ[9]) + 2 * r[10] * fabs(pQ[10]) + 2 * r[23] * fabs(pQ[23]) + 2 * r[25] * fabs(pQ[25]) + 2 * r[28] * fabs(pQ[28]));
		}
		
	
		/// add the adjust flow to each pQ pipes
	

	}
}
__global__ void Add_Q(float* pQ,float*dR ,int n)
{
	pQ[0] = pQ[0] + dR[0] + dR[9];
	pQ[1] = pQ[1] + dR[0] - dR[2];
	pQ[2] = pQ[2] + dR[2] - dR[4];
	pQ[3] = pQ[3] + dR[1] + dR[9];
	pQ[4] = pQ[4] - dR[1] + dR[2];
	pQ[5] = pQ[5] + dR[3] - dR[4];
	pQ[6] = pQ[6] + dR[0] - dR[1];
	pQ[7] = pQ[7] - dR[1] + dR[4];
	pQ[8] = pQ[8] - dR[4] + dR[5];
	pQ[9] = pQ[9] - dR[0] + dR[10];
	pQ[10] = pQ[10] - dR[3] + dR[10];
	pQ[11] = pQ[11] + dR[4] + dR[9];
	pQ[12] = pQ[12] - dR[0] + dR[3];
	pQ[13] = pQ[13] - dR[2] + dR[3];
	pQ[14] = pQ[14] + dR[3] - dR[5];
	pQ[15] = pQ[15] + dR[3] - dR[6];
	pQ[16] = pQ[16] + dR[6] - dR[5];
	pQ[17] = pQ[17] - dR[3] + dR[7];
	pQ[18] = pQ[18] - dR[6] + dR[7];
	pQ[19] = pQ[19] + dR[6] - dR[8];
	pQ[20] = pQ[20] - dR[5] + dR[8];
	pQ[21] = pQ[21] - dR[5];
	pQ[22] = pQ[22] + dR[8];
	pQ[23] = pQ[23] - dR[7] + dR[10];
	pQ[24] = pQ[24] + dR[7] - dR[8];
	pQ[25] = pQ[25] - dR[8] + dR[10];
	pQ[26] = pQ[26] + dR[9] + dR[10];
	pQ[27] = pQ[27] + dR[9];
	pQ[28] = pQ[28] + dR[10];


	
}