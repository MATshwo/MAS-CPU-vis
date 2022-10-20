#pragma once
#include<iostream>
#include"SeSchwarzPreconditioner.h"
#include<emmintrin.h>
#include "edge_face.h"
#include "read_sysAb.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> //for matrices
#include <glm/gtc/type_ptr.hpp>
#include "lvec.h"
#include<windows.h>
#include "main.h"
#pragma comment(lib,"winmm.lib")

//#include "main.cpp"
namespace my
{
	
	//SE::SeSchwarzPreconditioner MAS(float *dV,LargeVector<glm::vec3> &z0); //dV参数可以去掉，提供了一种传递参数的选择,该部分用于返回一个预处理的sh
	//void MAS(LargeVector<glm::mat3> sys_A, SE::SeSchwarzPreconditioner &sh);
	void MAS(LargeVector<glm::mat3> sys_A);
	//执行矩阵-向量乘法，结果存储在z0中
	//int Cmu(SE::SeSchwarzPreconditioner sh, LargeVector<glm::vec3>& z0,LargeVector<glm::vec3> sys_b);
}