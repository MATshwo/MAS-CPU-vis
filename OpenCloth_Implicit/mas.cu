#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <time.h>

//int main(void)
//{
//	// init data
//	int num = 10;
//	int a[10], b[10], c[10];
//	int* a_gpu, * b_gpu, * c_gpu;
//
//	for (int i = 0; i < num; i++)
//	{
//		a[i] = i;
//		b[i] = i * i;
//	}
//
//
//
//	cudaMalloc((void**)&a_gpu, num * sizeof(int));
//	cudaMalloc((void**)&b_gpu, num * sizeof(int));
//	cudaMalloc((void**)&c_gpu, num * sizeof(int));
//
//	// copy data
//	cudaMemcpy(a_gpu, a, num * sizeof(int), cudaMemcpyHostToDevice);
//	cudaMemcpy(b_gpu, b, num * sizeof(int), cudaMemcpyHostToDevice);
//
//	addKernel(a_gpu, b_gpu, c_gpu, num);
//
//	//printf("%d + %d = %d\n", a_gpu[0], b_gpu[0], c_gpu[0]);
//
//	// get data
//	cudaMemcpy(c, c_gpu, num * sizeof(int), cudaMemcpyDeviceToHost);
//
//	// visualization
//	for (int i = 0; i < num; i++)
//	{
//		printf("%d + %d = %d\n", a[i], b[i], c[i]);
//	}
//
//	return 0;
//}
//int main()
//{
//	getThreadNum();
//}