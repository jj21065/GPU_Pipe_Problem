
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <malloc.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<windows.h>
#define N 11
#define NN 29
#define TPB 29
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
int main()
{
	
	
	Allocate_Memory(N);
	cudaMemset(dR, 0.0,N*sizeof(float));

	float rs0 = 1;
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

	int n = 2;
	int iter_no = 5000;
	int i;

	printf("inital Q[0] = %g\n", Q[0]);

	LARGE_INTEGER startTime, endTime, fre;
	double times;
	QueryPerformanceFrequency(&fre); //取得CPU頻率
	QueryPerformanceCounter(&startTime); //取得開機到現在經過幾個CPU Cycle

	CopyMemToDevice(NN);
	for (i = 0; i < iter_no; i++)
	{
		Compute_Q<<<BPG,TPB>>>(d_Q,dR, N);
		Add_Q << <1, 1 >> >(d_Q, dR, N);
		
	}
	CopyMemToHost(NN);

	QueryPerformanceCounter(&endTime); //取得開機到程式執行完成經過幾個CPU Cycle
	times = ((double)endTime.QuadPart - (double)startTime.QuadPart) / fre.QuadPart;
	

	printf("time consume : %f\n", times);
	printf("final Q[0] = %g\n", Q[0]);
	printf("final Q9-EA = %g\n", Q[27]);
	printf("final Qh-EB = %g\n", Q[28]);

	system("pause");
	Free_Memory();
	return 0;
}

__global__ void Compute_Q(float* pQ,float*dR,int n )
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
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

	//float R[11] = { 0 };

	if (i < n)
	{
	
		// calculate adjust flow pQ
		if (i == 0)
		{
			if (2 * (r01*abs(pQ[1]) + r12*abs(pQ[3]) + r24*abs(pQ[6]) + r34*abs(pQ[5]) + r03*abs(pQ[4])) > 0)
				dR[0] = -(r01*pQ[1] * abs(pQ[1]) + r12*pQ[3] * abs(pQ[3]) + r24*pQ[6] * abs(pQ[6]) - r34*pQ[5] * abs(pQ[5]) - r03*pQ[4] * abs(pQ[4])) / (2 * (r01*abs(pQ[1]) + r12*abs(pQ[3]) + r24*abs(pQ[6]) + r34*abs(pQ[5]) + r03*abs(pQ[4])));
		}
		else if (i == 1)
		{
			if (2 * (r16*abs(pQ[2]) + r56*abs(pQ[8]) + r25*abs(pQ[7]) + r12*abs(pQ[3])) > 0)
				dR[1] = -(r16*pQ[2] * abs(pQ[2]) - r56*pQ[8] * abs(pQ[8]) - r25*pQ[7] * abs(pQ[7]) - r12*pQ[3] * abs(pQ[3])) / (2 * (r16*abs(pQ[2]) + r56*abs(pQ[8]) + r25*abs(pQ[7]) + r12*abs(pQ[3])));
		}
		else if (i == 2){
			if (2 * (r34*abs(pQ[5]) + r47*abs(pQ[10]) + r78*abs(pQ[12]) + r8a*abs(pQ[15]) + rab*abs(pQ[19]) + reb*abs(pQ[16]) + r3e*abs(pQ[9])) > 0)
				dR[2] = -(r34*pQ[5] * abs(pQ[5]) + r47*pQ[10] * abs(pQ[10]) + r78*pQ[12] * abs(pQ[12]) + r8a*pQ[15] * abs(pQ[15]) + rab*pQ[19] * abs(pQ[19]) - reb*pQ[16] * abs(pQ[16]) - r3e*pQ[9] * abs(pQ[9])) / (2 * (r34*abs(pQ[5]) + r47*abs(pQ[10]) + r78*abs(pQ[12]) + r8a*abs(pQ[15]) + rab*abs(pQ[19]) + reb*abs(pQ[16]) + r3e*abs(pQ[9])));
		}
		else if (i == 3){
			if (2 * (r56*abs(pQ[8]) + r69*abs(pQ[14]) + r89*abs(pQ[13]) + r78*abs(pQ[12]) + r57*abs(pQ[11])) > 0)
				dR[3] = -(r56*pQ[8] * abs(pQ[8]) + r69*pQ[14] * abs(pQ[14]) - r89*pQ[13] * abs(pQ[13]) - r78*pQ[12] * abs(pQ[12]) - r57*pQ[11] * abs(pQ[11])) / (2 * (r56*abs(pQ[8]) + r69*abs(pQ[14]) + r89*abs(pQ[13]) + r78*abs(pQ[12]) + r57*abs(pQ[11])));
		}
		else if (i == 4){
			if (2 * (r89*abs(pQ[13]) + r9g*abs(pQ[21]) + rdg*abs(pQ[17]) + rad*abs(pQ[23]) + r8a*abs(pQ[15])) > 0)
				dR[4] = -(r89*pQ[13] * abs(pQ[13]) - r9g*pQ[21] * abs(pQ[21]) - rdg*pQ[17] * abs(pQ[17]) - rad*pQ[23] * abs(pQ[23]) - r8a*pQ[15] * abs(pQ[15])) / (2 * (r89*abs(pQ[13]) + r9g*abs(pQ[21]) + rdg*abs(pQ[17]) + rad*abs(pQ[23]) + r8a*abs(pQ[15])));
		}
		else if (i == 5){
			if (2 * (rad*abs(pQ[23]) + rdc*abs(pQ[24]) + rbc*abs(pQ[20]) + rab*abs(pQ[19])) > 0)
				dR[5] = -(rad*pQ[23] * abs(pQ[23]) + rdc*pQ[24] * abs(pQ[24]) - rbc*pQ[20] * abs(pQ[20]) - rab*pQ[19] * abs(pQ[19])) / (2 * (rad*abs(pQ[23]) + rdc*abs(pQ[24]) + rbc*abs(pQ[20]) + rab*abs(pQ[19])));
		}
		else if (i == 6){
			if (2 * (reb*abs(pQ[16]) + rbc*abs(pQ[20]) + rcf*abs(pQ[22]) + ref*abs(pQ[18])) > 0)
				dR[6] = -(reb*pQ[16] * abs(pQ[16]) + rbc*pQ[20] * abs(pQ[20]) + rcf*pQ[22] * abs(pQ[22]) - ref*pQ[18] * abs(pQ[18])) / (2 * (reb*abs(pQ[16]) + rbc*abs(pQ[20]) + rcf*abs(pQ[22]) + ref*abs(pQ[18])));
		}
		else if (i == 7){
			if (2 * (rdg*abs(pQ[17]) + rgh*abs(pQ[25]) + rfh*abs(pQ[26]) + rcf*abs(pQ[22]) + rdc*abs(pQ[24])) > 0)
				dR[7] = -(rdg*pQ[17] * abs(pQ[17]) + rgh*pQ[25] * abs(pQ[25]) - rfh*pQ[26] * abs(pQ[26]) - rcf*pQ[22] * abs(pQ[22]) - rdc*pQ[24] * abs(pQ[24])) / (2 * (rdg*abs(pQ[17]) + rgh*abs(pQ[25]) + rfh*abs(pQ[26]) + rcf*abs(pQ[22]) + rdc*abs(pQ[24])));
		}
		else if (i == 8)
		{
			if (2 * (r25*abs(pQ[7]) + r57*abs(pQ[11]) + r47*abs(pQ[10]) + r24*abs(pQ[6])) > 0)
				dR[8] = -(r25*pQ[7] * abs(pQ[7]) + r57*pQ[11] * abs(pQ[11]) - r47*pQ[10] * abs(pQ[10]) - r24*pQ[6] * abs(pQ[6])) / (2 * (r25*abs(pQ[7]) + r57*abs(pQ[11]) + r47*abs(pQ[10]) + r24*abs(pQ[6])));
		}
		else if (i == 9){
			if (2 * (rs0*abs(pQ[0]) + r01*abs(pQ[1]) + r16*abs(pQ[2]) + r69*abs(pQ[14]) + r9EA*abs(pQ[27])) > 0)
				dR[9] = -(rs0*pQ[0] * abs(pQ[0]) + r01*pQ[1] * abs(pQ[1]) + r16*pQ[2] * abs(pQ[2]) + r69*pQ[14] * abs(pQ[14]) + r9EA*pQ[27] * abs(pQ[27]) - 10) / (2 * (rs0*abs(pQ[0]) + r01*abs(pQ[1]) + r16*abs(pQ[2]) + r69*abs(pQ[14]) + r9EA*abs(pQ[27])));
		}
		else if (i == 10){
			if (2 * (rs0*abs(pQ[0]) + r03*abs(pQ[4]) + r3e*abs(pQ[9]) + ref*abs(pQ[18]) + rfh*abs(pQ[26]) + rhEB*abs(pQ[28])) > 0)
				dR[10] = -(rs0*pQ[0] * abs(pQ[0]) + r03*pQ[4] * abs(pQ[4]) + r3e*pQ[9] * abs(pQ[9]) + ref*pQ[18] * abs(pQ[18]) + rfh*pQ[26] * abs(pQ[26]) + rhEB*pQ[28] * abs(pQ[28]) - 10) / (2 * (rs0*abs(pQ[0]) + r03*abs(pQ[4]) + r3e*abs(pQ[9]) + ref*abs(pQ[18]) + rfh*abs(pQ[26]) + rhEB*abs(pQ[28])));
		}
		
	
		/// add the adjust flow to each pQ pipes
	

	}
}
__global__ void Add_Q(float* pQ,float*dR ,int n)
{
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	pQ[1] = pQ[1] + dR[0] + dR[9];
	pQ[3] = pQ[3] + dR[0] - dR[1];
	pQ[6] = pQ[6] + dR[0] - dR[8];
	pQ[5] = pQ[5] - dR[0] + dR[2];
	pQ[4] = pQ[4] - dR[0] + dR[10];
	pQ[2] = pQ[2] + dR[1] + dR[9];
	pQ[8] = pQ[8] - dR[1] + dR[3];
	pQ[7] = pQ[7] - dR[1] + dR[8];
	pQ[10] = pQ[10] - dR[8] + dR[2];
	pQ[12] = pQ[12] + dR[2] - dR[3];
	pQ[15] = pQ[15] + dR[2] - dR[4];

	pQ[11] = pQ[11] - dR[3] + dR[8];
	pQ[13] = pQ[13] - dR[3] + dR[4];
	pQ[14] = pQ[14] + dR[3] + dR[9];
	pQ[23] = pQ[23] - dR[4] + dR[5];
	pQ[21] = pQ[21] - dR[4];
	pQ[17] = pQ[17] - dR[4] + dR[7];
	pQ[24] = pQ[24] + dR[5] - dR[7];
	pQ[25] = pQ[25] + dR[7];
	pQ[26] = pQ[26] - dR[7] + dR[10];
	pQ[28] = pQ[28] + dR[10];
	pQ[27] = pQ[27] + dR[9];
	pQ[19] = pQ[19] + dR[2] - dR[5];
	pQ[16] = pQ[16] - dR[2] + dR[6];
	pQ[9] = pQ[9] - dR[2] + dR[10];
	pQ[18] = pQ[18] - dR[6] + dR[10];
	pQ[20] = pQ[20] - dR[5] + dR[6];
	pQ[22] = pQ[22] + dR[6] - dR[7];
	pQ[0] = pQ[0] + dR[10] + dR[9];
}