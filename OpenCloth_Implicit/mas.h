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
	
	//SE::SeSchwarzPreconditioner MAS(float *dV,LargeVector<glm::vec3> &z0); //dV��������ȥ�����ṩ��һ�ִ��ݲ�����ѡ��,�ò������ڷ���һ��Ԥ�����sh
	//void MAS(LargeVector<glm::mat3> sys_A, SE::SeSchwarzPreconditioner &sh);
	void MAS(LargeVector<glm::mat3> sys_A);
	//ִ�о���-�����˷�������洢��z0��
	//int Cmu(SE::SeSchwarzPreconditioner sh, LargeVector<glm::vec3>& z0,LargeVector<glm::vec3> sys_b);
}