
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

float mypow(float, int);
float myabs(float a)
{
	if (a < 0)
		return -a;
	else
		return a;
}
int main()
{
	// pipe parameter
	float r12 = 5;
	float r13 = 1;
	float r23 = 10;
	float r24 = 1;
	float r34 = 5;
	// initial conditions
	float Q12 = 0;
	float Q13 = 10;
	float Q23 = -3;
	float Q24 = Q12 - Q23;
	float Q34 = Q13 + Q23;

	int n = 2;
	int iter_no = 500;
	int i;
	//r23 = (r24*Q24*Q24-r34*Q34*Q34)/Q23/Q23;
	//r23 = (r12*Q12*Q12-r13*Q13*Q13)/Q23/Q23;
	int c12 = 1, c13 = 1, c23 = 1, c24 = 1, c34 = 1;

	for (i = 0; i < iter_no; i++)
	{
		//Control the direction sign
		c12 = (Q12 >= 0) ? 1 : -1;
		c13 = (Q13 >= 0) ? 1 : -1;
		c23 = (Q23 >= 0) ? 1 : -1;
		c24 = (Q24 >= 0) ? 1 : -1;
		c34 = (Q34 >= 0) ? 1 : -1;
		float dQ_2 = -(r24*c24*Q24*Q24 - c34*r34*Q34*Q34 - c23*r23*Q23*Q23) / (2 * (r24*myabs(Q24) + r34*myabs(Q34) + r23*myabs(Q23)));
		Q12 = Q12 + dQ_2;
		Q13 = Q13 - dQ_2;
		Q24 = Q24 + dQ_2;
		Q34 = Q34 - dQ_2;
		//calculate the new r23 from the last Q
		r23 = myabs((-r12*c12*Q12*Q12 + c13*r13*Q13*Q13) / (Q23*Q23*c23));
	}


	printf("\n\n");
	printf("Q12 = %g \n", Q12);
	printf("Q13 = %g \n", Q13);
	printf("Q23 = %g \n", Q23);
	printf("Q24 = %g \n", Q24);
	printf("Q34 = %g \n\n", Q34);
	printf("r23 = %g \n", r23);
	system("pause");
	return 0;
}
float mypow(float value, int n)
{
	int i;
	if (n == 0)
		return 1;
	for (i = 0; i < n - 1; i++)
	{
		value = value*value;
	}
	return value;
}
