/*
Copyright (c) 2011, Movania Muhammad Mobeen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

//A simple cloth using implicit integration based on the model by Baraff
//& Witkin in the paper "Large steps in cloth simulation" and the SIGGRAPH course notes
//"Realtime Physics" http://www.matthiasmueller.info/realtimephysics/coursenotes.pdf using
//GLUT,GLEW and GLM libraries. This code is intended for beginners so that they may
//understand what is required to implement implicit integration based cloth simulation.
//
//This code is under BSD license. If you make some improvements,
//or are using this in your research, do let me know and I would appreciate
//if you acknowledge this in your code.
//
//Controls:
//left click on any empty region to rotate, middle click to zoom
//left click and drag any point to drag it.
//
//Author: Movania Muhammad Mobeen
//        School of Computer Engineering,
//        Nanyang Technological University,
//        Singapore.
//Email : mova0002@e.ntu.edu.sg


/*MAS+PCG复现+优化后的版本*/

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> //for matrices
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
#include "mas.h"
#include "lvec.h"
//#include<windows.h>
//#pragma comment(lib,"winmm.lib")
#pragma comment(lib, "glew32.lib")

// 全局变量里加上这几个
// 鼠标状态
bool isLeftButtonDown = false;
bool isRighttButtonDown = false;
bool isMidtButtonDown = false;

SE::SeSchwarzPreconditioner sh;
int diagflag = 0; //用于指定MAS何时重新计算
int flag = 0; //统计MAS预计算的次数

char A_path[] = "./Hessian_matrix1.txt";
char b_path[] = "./Residual1.txt";
using namespace std;
//窗口大小
const int width = 1024, height = 1024;

//顶点个数
int numX = 300, numY= 300; //
const int total_points = (numX+1)*(numY+1);
float fullsize = 4.0f;
float halfsize = fullsize/2.0f;

//仿真时间和速率
float timeStep =  10/60.0f; //时间步长
float currentTime = 0;
double accumulator = timeStep;
int selected_index = -1;

struct Spring {
	int p1, p2; //端点
	float rest_length; //静息长度
	float Ks, Kd; //刚度系数和阻尼系数
	int type; 
};

vector<GLushort> indices; //三角形面片索引：每三个构成一个三角形
vector<Spring> springs;

vector<glm::vec3> X;
vector<glm::vec3> V;
vector<glm::vec3> F;
vector<glm::vec3> dc_dp; //  dc/dp c关于位置的导数
vector<glm::mat3> df_dx; //  df/dp 合外力关于位置的导数
vector<glm::mat3> df_dv; //  df/dv 外力关于速度的导数
vector<float> C; //for implicit integration 浮点型的向量
vector<float> C_Dot; //for implicit integration 文中的一种写法：具有特定含义
vector<glm::vec3> deltaP2; 

char info[MAX_PATH]={0};

int oldX=0, oldY=0;
float rX=15, rY=0; //摄像机的位置/观察角度
int state =1 ;
float dist=-30; //摄像机距离布料的距离(负值)，正值图像会消失
const int GRID_SIZE=10;

const int STRUCTURAL_SPRING = 0;
const int SHEAR_SPRING = 1;
const int BEND_SPRING = 2;

int spring_count=0;

const float DEFAULT_DAMPING =  -0.125f;
float	KsStruct = 1.75f,KdStruct = -0.25f;
float	KsShear = 1.75f,KdShear = -0.25f;
float	KsBend = 1.95f,KdBend = -0.25f;
glm::vec3 gravity=glm::vec3(0.0f,-0.00981f,0.0f);
float mass =1.0f;

GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;
float startTime =0, fps=0;
int totalFrames=0;
int totalTime =0;
	
const float EPS = 0.001f;
const float EPS2 = EPS*EPS;
const int i_max = 10;

glm::mat4 ellipsoid, inverse_ellipsoid;
int iStacks = 30;
int iSlices = 30;
float fRadius = 1;

// Resolve constraint in object space
glm::vec3 center = glm::vec3(0,0,0); //object space center of ellipsoid
float radius = 1;					 //object space radius of ellipsoid

void StepPhysics(float dt );
void SolveConjugateGradientPreconditioned(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b, LargeVector<glm::vec3> P, LargeVector<glm::vec3> P_inv) {
	//预处理的CG求解器 
	LARGE_INTEGER t0, tt, tc;
	QueryPerformanceFrequency(&tc);//tc是频率，本质上是计算系统打点的次数除以频率
	QueryPerformanceCounter(&t0);//起始时间
	int canPrintInfo = 0;

	float i = 0;
	LargeVector<glm::vec3> r = (b - A * x);
	LargeVector<glm::vec3> d = P_inv * r;
	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta = 0;
	float delta_old = 0;
	//float delta_new = dot(r,P*r); //按照cg算法流程图这里是p-1才对把
	float delta_new = dot(r, P_inv * r);
	float delta0 = delta_new;
	while (i<i_max && delta_new> EPS2 * delta0) {
		//clock_t time_stt = clock();
		q = A * d;
		alpha = delta_new / dot(d, q);
		x = x + alpha * d;
		r = r - alpha * q;
		delta_old = delta_new;
		delta_new = dot(r, P_inv * r); //一步到位？
		//delta_new = dot(r, r);
		beta = delta_new / delta_old;
		d = P_inv * r + beta * d;
		//d = r + beta * d;
		i++;

	}
	if (canPrintInfo)//可以传个bool变量进去，每隔100帧再输出
	{
		QueryPerformanceCounter(&tt);//当前时间
		std::cout << "Time Cost in JacobiPCG " << (double)((tt.QuadPart - t0.QuadPart) * 1000.0 / tc.QuadPart) << "ms" << std::endl;
		//std::cout << "PreparePreconditioner run success!" << std::endl;
	}
	
}

//void SolveConjugateGradientMAS(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b, SE::SeSchwarzPreconditioner sh) {
void SolveConjugateGradientMAS(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b) {
	
	LARGE_INTEGER t0, tt, tc;
	QueryPerformanceFrequency(&tc);//tc是频率，本质上是计算系统打点的次数除以频率
	QueryPerformanceCounter(&t0);//起始时间
	int canPrintInfo = 0;

	float i = 0;
	LargeVector<glm::vec3> r = (b - A * x);
	LargeVector<glm::vec3> z;
	z.resize(total_points);

	SE::SeVec3fSimd* z1 = new SE::SeVec3fSimd[sh.m_numVerts];
	SE::SeVec3fSimd* residual = new SE::SeVec3fSimd[sh.m_numVerts];
	for (int a = 0; a < sh.m_numVerts; a++)
	{
		//float temp[3] = { 0.1,0.2,0.3 };
		SE::SeVec3fSimd nows(r[a][0], r[a][1], r[a][2]); 
		residual[a] = nows;
	}
	const int dim = 1;
	sh.Preconditioning(z1, residual, dim);

	for (int k = 0; k < sh.m_numVerts; k++)
	{
		z[k] = glm::vec3(z1[k].x, z1[k].y, z1[k].z);
	}
	delete z1;
	delete residual;
	LargeVector<glm::vec3> d = z; //p0
	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta = 0;
	float delta_old = 0;
	float delta_new = dot(r,z); //delta_new = dot(r,z); 初始分子
	float delta0 = delta_new;
	while (i<i_max && delta_new> EPS2 * delta0) {
		q = A * d;
		alpha = delta_new /dot(d,q);
		x = x + alpha * d;
		r = r - alpha * q; 

		SE::SeVec3fSimd* z1 = new SE::SeVec3fSimd[sh.m_numVerts]; 
		SE::SeVec3fSimd* residual = new SE::SeVec3fSimd[sh.m_numVerts];
		for (int a = 0; a < sh.m_numVerts; a++)
		{
			SE::SeVec3fSimd nows(r[a][0], r[a][1], r[a][2]);
			residual[a] = nows;
		}
		const int dim = 1;
		sh.Preconditioning(z1, residual, dim);
		////std::cout << "Preconditioning run success!" << std::endl;
		//// 40次 42ms
		for (int k = 0; k < sh.m_numVerts; k++)
		{
			z[k] = glm::vec3(z1[k].x, z1[k].y, z1[k].z);
		}
		delta_old = delta_new;
		delta_new = dot(z, r);
		beta = delta_new / delta_old; //求出系数β
		d = r + beta * d; //d = z + beta * d; 更新p
		i++;
		delete z1;
		delete residual;
		//std::cout << "Now times is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << std::endl;
	}
	
	if (canPrintInfo)//可以传个bool变量进去，每隔100帧再输出
	{
		QueryPerformanceCounter(&tt);//当前时间
		std::cout << "Time Cost in MPCG " << (double)((tt.QuadPart - t0.QuadPart) * 1000.0 / tc.QuadPart) << "ms" << std::endl;
		//std::cout << "PreparePreconditioner run success!" << std::endl;
	}
	
}

LargeVector<glm::vec3> dV;
LargeVector<glm::mat3> A; 
glm::mat3 M = glm::mat3(1.0f); //3*3的单位矩阵
LargeVector<glm::vec3> b;
LargeVector<glm::vec3> P_;
LargeVector<glm::vec3> P_inv;
vector<float> inv_len;

void AddSpring(int a, int b, float ks, float kd, int type) {
	Spring spring;
	spring.p1=a;
	spring.p2=b;
	spring.Ks=ks;
	spring.Kd=kd;
	spring.type = type;
	glm::vec3 deltaP = X[a]-X[b];
	spring.rest_length = sqrt(glm::dot(deltaP, deltaP));
	springs.push_back(spring);
}
void OnMouseDown(int button, int s, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		if (s == GLUT_DOWN)
		{
			isLeftButtonDown = true;
			oldX = x;
			oldY = y;
			int window_y = (height - y);
			float norm_y = float(window_y) / float(height / 2.0);
			int window_x = x;
			float norm_x = float(window_x) / float(width / 2.0);

			float winZ = 0;
			glReadPixels(x, height - y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
			if (winZ == 1)
				winZ = 0;
			double objX = 0, objY = 0, objZ = 0;
			gluUnProject(window_x, window_y, winZ, MV, P, viewport, &objX, &objY, &objZ);
			glm::vec3 pt(objX, objY, objZ);
			size_t i = 0;
			for (i = 0; i < total_points; i++) {
				if (glm::distance(X[i], pt) < 0.1) {
					selected_index = i;
					//printf("Intersected at %d\n",i);
					break;
				}
			}
		}
		if (s == GLUT_UP)
		{
			isLeftButtonDown = false;

			selected_index = -1;
			glutSetCursor(GLUT_CURSOR_INHERIT);
		}
	}

	else if (button == GLUT_MIDDLE_BUTTON)
	{
		if (s == GLUT_DOWN) isMidtButtonDown = true;
		if (s == GLUT_UP) isMidtButtonDown = false;
	}

	else if (button == GLUT_RIGHT_BUTTON)
	{
		if (s == GLUT_DOWN)
		{
			isRighttButtonDown = true;
			oldX = x;
			oldY = y;
		}

		if (s == GLUT_UP) isRighttButtonDown = false;
	}
	// 加了个滚轮事件
	else if (button == 3 && s == GLUT_DOWN)
	{
		dist *= (1 + (-10) / 60.0f);
	}
	else if (button == 4 && s == GLUT_DOWN)
	{
		dist *= (1 + (10) / 60.0f);
	}
}


void OnMouseMove(int x, int y)
{
	if (selected_index == -1) {//如果未选中点
		if (isRighttButtonDown)
		{
			rY += (x - oldX) / 5.0f;
			rX += (y - oldY) / 5.0f;
		}
	}
	else {//拖动选中点
		float delta = 1500 / abs(dist);
		float valX = (x - oldX) / delta;
		float valY = (oldY - y) / delta;
		if (abs(valX) > abs(valY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else
			glutSetCursor(GLUT_CURSOR_UP_DOWN);

		V[selected_index] = glm::vec3(0);
		X[selected_index].x += Right[0] * valX;
		float newValue = X[selected_index].y + Up[1] * valY;
		if (newValue > 0)
			X[selected_index].y = newValue;
		X[selected_index].z += Right[2] * valX + Up[2] * valY;
	}
	oldX = x;
	oldY = y;
	glutPostRedisplay();
}

void DrawGrid()
{
	//地面网格绘制：绘制地面的网格线条
	glBegin(GL_LINES);
	glColor3f(0.1f, 0.5f, 0.5f);
	for(int i=-GRID_SIZE;i<=GRID_SIZE;i++)
	{
		glVertex3f((float)i,0,(float)-GRID_SIZE);
		glVertex3f((float)i,0,(float)GRID_SIZE);

		glVertex3f((float)-GRID_SIZE,0,(float)i);
		glVertex3f((float)GRID_SIZE,0,(float)i);
	}
	glEnd();
}

void InitGL() 
{
	//initialize

	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	currentTime = startTime;

	// get ticks per second
    QueryPerformanceFrequency(&frequency);

    // start timer
    QueryPerformanceCounter(&t1);

	glEnable(GL_DEPTH_TEST);
	int i=0, j=0, count=0;
	int l1=0, l2=0;
	float ypos = 7.0f;
	int v = numY+1;
	int u = numX+1;

	indices.resize( numX*numY*2*3); //total-node = numx+1)*(numy+1
	X.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);

	A.resize(total_points); 
	b.resize(total_points);
	dV.resize(total_points);
	P_.resize(total_points);
	P_inv.resize(total_points);


	//fill in X
	for( j=0;j<=numY;j++) {
		for( i=0;i<=numX;i++) {
			X[count++] = glm::vec3( ((float(i)/(u-1)) *2-1)* halfsize, fullsize+1, ((float(j)/(v-1) )* fullsize));
		}
	}

	//fill in V 初始=0
	memset(&(V[0].x),0,total_points*sizeof(glm::vec3));

	//fill in indices 
	GLushort* id=&indices[0];
	for (i = 0; i < numY; i++) {
		for (j = 0; j < numX; j++) {
			int i0 = i * (numX+1) + j;
			int i1 = i0 + 1;
			int i2 = i0 + (numX+1);
			int i3 = i2 + 1;
			if ((j+i)%2) {
				*id++ = i0; *id++ = i2; *id++ = i1;
				*id++ = i1; *id++ = i2; *id++ = i3;
			} else {
				*id++ = i0; *id++ = i2; *id++ = i3;
				*id++ = i0; *id++ = i3; *id++ = i1;
			}
		}
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glPolygonMode(GL_BACK, GL_LINE);
	glPointSize(5);

	wglSwapIntervalEXT(0);

	//setup springs
	// Horizontal 水平的拉伸弹簧
	for (l1 = 0; l1 < v; l1++)	// v:垂直方向的顶点个数
		for (l2 = 0; l2 < (u - 1); l2++) {
			AddSpring((l1 * u) + l2,(l1 * u) + l2 + 1,KsStruct,KdStruct,STRUCTURAL_SPRING);
		}

	// Vertical 垂直方向的拉伸弹簧  u:水平方向的顶点个数
	for (l1 = 0; l1 < (u); l1++)
		for (l2 = 0; l2 < (v - 1); l2++) {
			AddSpring((l2 * u) + l1,((l2 + 1) * u) + l1,KsStruct,KdStruct,STRUCTURAL_SPRING);
		}


	// Shearing Springs 每个三角形两条对角线
	for (l1 = 0; l1 < (v - 1); l1++)
		for (l2 = 0; l2 < (u - 1); l2++) {
			AddSpring((l1 * u) + l2,((l1 + 1) * u) + l2 + 1,KsShear,KdShear,SHEAR_SPRING);
			AddSpring(((l1 + 1) * u) + l2,(l1 * u) + l2 + 1,KsShear,KdShear,SHEAR_SPRING);
		}

	// Bend Springs
	for (l1 = 0; l1 < (v); l1++) {
		for (l2 = 0; l2 < (u - 2); l2++) {
			AddSpring((l1 * u) + l2,(l1 * u) + l2 + 2,KsBend,KdBend,BEND_SPRING);
		}
		AddSpring((l1 * u) + (u - 3),(l1 * u) + (u - 1),KsBend,KdBend,BEND_SPRING); //这里的弹簧与循环上边最后一条弹簧重叠了？
	}
	for (l1 = 0; l1 < (u); l1++) {
		for (l2 = 0; l2 < (v - 2); l2++) {
			AddSpring((l2 * u) + l1,((l2 + 2) * u) + l1,KsBend,KdBend,BEND_SPRING);
		}
		AddSpring(((v - 3) * u) + l1,((v - 1) * u) + l1,KsBend,KdBend,BEND_SPRING);//为什么要重复
	}

	//初始化参数：用0填充
	int total_springs = springs.size();
	C.resize(total_springs ); //条件矩阵C(x)
	inv_len.resize(total_springs );
	C_Dot.resize(total_springs );
	dc_dp.resize(total_springs );
	deltaP2.resize(total_springs );
	df_dx.resize(total_springs );
	df_dv.resize(total_springs );
	memset(&(C[0]),0,total_springs*sizeof(float));
	memset(&(C_Dot[0]),0,total_springs*sizeof(float));
	memset(&(deltaP2[0].x),0,total_springs*sizeof(glm::vec3));

 	memset(&(P_[0].x),0,total_points*sizeof(glm::vec3));
	memset(&(P_inv[0].x),0,total_points*sizeof(glm::vec3));
	
	//create a basic ellipsoid object
	ellipsoid = glm::translate(glm::mat4(1),glm::vec3(0,2,0));
	ellipsoid = glm::rotate(ellipsoid, 45.0f ,glm::vec3(1,0,0));
	ellipsoid = glm::scale(ellipsoid, glm::vec3(fRadius,fRadius,fRadius/2));
	inverse_ellipsoid = glm::inverse(ellipsoid);
}

void OnReshape(int nw, int nh) {

	glViewport(0,0,nw, nh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(40, (GLfloat)nw / (GLfloat)nh, 1.f, 100.0f); //第一个参数是观察距离。后两个是长度宽度，最后两个参数猜测指定了一个范围

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_PROJECTION_MATRIX, P);

	glMatrixMode(GL_MODELVIEW);
}

void OnRender() {

	//时间和帧设置
	size_t i=0;
	float newTime = (float) glutGet(GLUT_ELAPSED_TIME);
	frameTime = newTime-currentTime; //
	currentTime = newTime;
	//accumulator += frameTime;

	//Using high res. counter
    QueryPerformanceCounter(&t2);
	 // compute and print the elapsed time in millisec
    frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
	t1=t2;
	accumulator += frameTimeQP;

	++totalFrames;
	if((newTime-startTime)>1000)
	{
		float elapsedTime = (newTime-startTime);
		fps = (totalFrames/ elapsedTime)*1000 ;
		startTime = newTime;
		totalFrames=0;
	}
	//在窗口显示模拟过程的一些参数：帧速率等
	//sprint 将格式化的数据输入到某个字符串
	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f", fps, frameTime, frameTimeQP);
	std::cout << "FPS:" << fps << std::endl;
	
	glutSetWindowTitle(info); //在显示窗口的标题处显示定义的字符串
	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
	
	//设置照相机属性
	glLoadIdentity();
	glTranslatef(0,0,dist);
	glRotatef(rX,1,0,0); 
	glRotatef(rY,0,1,0); 

	glGetDoublev(GL_MODELVIEW_MATRIX, MV);
	viewDir.x = (float)-MV[2];
	viewDir.y = (float)-MV[6];
	viewDir.z = (float)-MV[10];
	Right = glm::cross(viewDir, Up);

	DrawGrid();

	//draw ellipsoid
	glColor3f(0,0.4,0.3);
	glPushMatrix(); 
	glMultMatrixf(glm::value_ptr(ellipsoid));
	glutWireSphere(fRadius, iSlices, iStacks);
	glPopMatrix(); 

	//draw polygons
	glColor3f(0.1f,0.6f,0.0f);
	glBegin(GL_TRIANGLES);

	for(i=0;i<indices.size();i+=3) {
		glm::vec3 p1 = X[indices[i]];
		glm::vec3 p2 = X[indices[i+1]];
		glm::vec3 p3 = X[indices[i+2]];
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
		glVertex3f(p3.x,p3.y,p3.z);
	}
	glEnd();

	//draw points
	glBegin(GL_POINTS);
	for(i=0;i<total_points;i++) {
		glm::vec3 p = X[i];
		int is = (i==selected_index);
		glColor3f(0.6f, 0.3f, 0.3f);
		//glColor3f((float)!is,(float)is,(float)is);
		glVertex3f(p.x,p.y,p.z);
	}
	glEnd();

	glutSwapBuffers();
}

void OnShutdown() {
	X.clear();
	V.clear();
	F.clear();
	indices.clear();
	springs.clear();
	dc_dp.clear();
	df_dx.clear();
	df_dv.clear();
	C.clear();
	C_Dot.clear();
	deltaP2.clear();

	dV.clear();
	A.clear();

	b.clear();
	P_.clear();
	P_inv.clear();
	inv_len.clear();
}

void ComputeForces() {
	size_t i=0;
	for (i = 0; i < total_points; i++) {
		F[i] = glm::vec3(0);
	
		//add gravity force
		if (i != 0 && i != (numX))
			F[i] += gravity;
		//add force due to damping of velocity 阻尼关于速度的函数
		F[i] += DEFAULT_DAMPING * V[i];
	}
	
	for (i = 0; i < springs.size(); i++) {
		glm::vec3 p1 = X[springs[i].p1];
		glm::vec3 p2 = X[springs[i].p2];
		glm::vec3 v1 = V[springs[i].p1];
		glm::vec3 v2 = V[springs[i].p2];
		glm::vec3 deltaP = p1 - p2;
		glm::vec3 deltaV = v1 - v2;
		float dist = glm::length(deltaP);
		inv_len[i] = 1.0f / dist;
		C[i] = dist - springs[i].rest_length;
		dc_dp[i] = deltaP / dist;
		C_Dot[i] = glm::dot(v1, -dc_dp[i]) + glm::dot(v2, dc_dp[i]);
		deltaP2[i] = glm::vec3(deltaP.x * deltaP.x,deltaP.y * deltaP.y, deltaP.z * deltaP.z);
	
		float leftTerm = -springs[i].Ks * (dist - springs[i].rest_length);
		float rightTerm = +springs[i].Kd * (glm::dot(deltaV, deltaP) / dist);
		glm::vec3 springForce = (leftTerm + rightTerm) * glm::normalize(deltaP);
	
		if (springs[i].p1 != 0 && springs[i].p1 != numX)
			F[springs[i].p1] += springForce;
		if (springs[i].p2 != 0 && springs[i].p2 != numX)
			F[springs[i].p2] -= springForce;
	}
	
}

void CalcForceDerivatives() {
	//clear the derivatives
	memset(&(df_dx[0]),0,total_points*sizeof(glm::mat3));
	memset(&(df_dv[0]),0,total_points*sizeof(glm::mat3));

	size_t i=0;

	glm::mat3 d2C_dp2[2][2]={glm::mat3(1.0f),glm::mat3(1.0f),glm::mat3(1.0f),glm::mat3(1.0f)}; //求二阶导


	for (i = 0; i < springs.size(); i++) {
		float c1 = C[i];
		d2C_dp2[0][0][0][0] = (-c1 * deltaP2[i].x + c1);
		d2C_dp2[0][0][1][1] = (-c1 * deltaP2[i].y + c1);
		d2C_dp2[0][0][2][2] = (-c1 * deltaP2[i].z + c1);

		d2C_dp2[0][1][0][0] = (c1 * deltaP2[i].x - c1);
		d2C_dp2[0][1][1][1] = (c1 * deltaP2[i].y - c1);
		d2C_dp2[0][1][2][2] = (c1 * deltaP2[i].z - c1);

		d2C_dp2[1][0] = d2C_dp2[0][1];
		d2C_dp2[1][1] = d2C_dp2[0][0];

		glm::mat3 dp1 = glm::outerProduct(dc_dp[i], dc_dp[i]); //外积：向量乘法得到一个矩阵
		glm::mat3 dp2 = glm::outerProduct(dc_dp[i], -dc_dp[i]);
		glm::mat3 dp3 = glm::outerProduct(-dc_dp[i],-dc_dp[i]);

		//df/dx只考虑弹簧和阻尼的影响
		df_dx[i] += -springs[i].Ks * (dp1 + (d2C_dp2[0][0] * C[i])) - springs[i].Kd * (d2C_dp2[0][0] * C_Dot[i]);
		df_dx[i] += -springs[i].Ks * (dp2 + (d2C_dp2[0][1] * C[i])) - springs[i].Kd * (d2C_dp2[0][1] * C_Dot[i]);
		df_dx[i] += -springs[i].Ks * (dp2 + (d2C_dp2[1][1] * C[i])) - springs[i].Kd * (d2C_dp2[1][1] * C_Dot[i]);

		//df/dv主要是阻尼的影响
		df_dv[i] += -springs[i].Kd * dp1;
		df_dv[i] += -springs[i].Kd * dp2;
		df_dv[i] += -springs[i].Kd * dp3;
	}
	
}


void IntegrateImplicit(float deltaTime) {
	

	LARGE_INTEGER t0, tt, tc;
	QueryPerformanceFrequency(&tc);//tc是频率，本质上是计算系统打点的次数除以频率
	QueryPerformanceCounter(&t0);//起始时间
	int canPrintInfo = 1;


	float h = deltaTime;
	CalcForceDerivatives(); //
    float y = 0.0;//correction term
	size_t i=0;

	
	for (i = 0; i < total_points; i++) {

		//系统矩阵的两种选择：1）牛顿法：能量函数的二阶导数；2）梯度下降法：单位矩阵
		A[i] = M - h * (df_dv[i] + h * df_dx[i]); 
		//A[i] = glm::mat3(1,0,0,0,1,0,0,0,1);
		b[i] = h * (F[i] + h * df_dx[i] * (V[i] + y));
		P_[i] = glm::vec3(A[i][0][0], A[i][1][1], A[i][2][2]);
		P_inv[i] = 1.0f / P_[i];
	}
	

	// //Jacobi+PCG
	//SolveConjugateGradientPreconditioned(A, dV, b, P_, P_inv); 


	
	////if (diagflag == 0)  { my::MAS(A, sh); diagflag = -1; } //错误代码不用管
	

	//隐式积分+MAS
	//my::MAS(A);   //每次执行---功能：根据顶点坐标等相关信息预计算MAS处理器，对应AllocatePrecoditioner()&PreparePreconditioner()部分
	//if (diagflag == 0)  { my::MAS(A); diagflag = -1; }//整个仿真阶段只执行一次MAS的预计算
	if (flag % 10 == 0)   { my::MAS(A); } //printf("____________________%d__________________\n", flag); } //每隔10个时间步更新一次MAS
	flag += 1;
	
	//printf("%d\n", sh.m_numVerts);
	//sh(预处理器)定义为全局变量更合适，省去了参数传递的耗时(sh作为参数传递耗时大概10-20ms)
	SolveConjugateGradientMAS(A,dV,b); //wanghuamin MAS+PCG

	for(i=0;i<total_points;i++) {
		V[i] += (dV[i]*deltaTime);
		X[i] += deltaTime*V[i];
		if(X[i].y <0) 
		{
			X[i].y = 0;
		}
	}

	if (canPrintInfo)//可以传个bool变量进去，每隔100帧再输出
	{
		QueryPerformanceCounter(&tt);//当前时间
		std::cout << "Total time cost per timeStep: " << (double)((tt.QuadPart - t0.QuadPart) * 1000.0 / tc.QuadPart) << "ms" << std::endl;
		//std::cout << "PreparePreconditioner run success!" << std::endl;
	}
	
}

void ApplyProvotDynamicInverse() {

	for (size_t i = 0; i < springs.size(); i++) {
		//check the current lengths of all springs
		glm::vec3 p1 = X[springs[i].p1];
		glm::vec3 p2 = X[springs[i].p2];
		glm::vec3 deltaP = p1 - p2;
		float dist = glm::length(deltaP);
		if (dist > (springs[i].rest_length * 1.01f))
		{

			dist -= (springs[i].rest_length * 1.01f);
			dist /= 2.0f;
			deltaP = glm::normalize(deltaP);
			deltaP *= dist;

			if (springs[i].p1 == 0 || springs[i].p1 == numX) {
				V[springs[i].p2] += deltaP;
			}
			else if (springs[i].p2 == 0 || springs[i].p2 == numX) {
				V[springs[i].p1] -= deltaP;
			}
			else {
				V[springs[i].p1] -= deltaP;
				V[springs[i].p2] += deltaP;
			}
		}
	}
}

void EllipsoidCollision() {
	for(size_t i=0;i<total_points;i++) {
		glm::vec4 X_0 = (inverse_ellipsoid*glm::vec4(X[i],1));
		glm::vec3 delta0 = glm::vec3(X_0.x, X_0.y, X_0.z) - center;
		float distance = glm::length(delta0);
		if (distance < 1.0f) {
			delta0 = (radius - distance) * delta0 / distance;

			// Transform the delta back to original space
			glm::vec3 delta;
			glm::vec3 transformInv;
			transformInv = glm::vec3(ellipsoid[0].x, ellipsoid[1].x, ellipsoid[2].x);
			transformInv /= glm::dot(transformInv, transformInv);
			delta.x = glm::dot(delta0, transformInv);
			transformInv = glm::vec3(ellipsoid[0].y, ellipsoid[1].y, ellipsoid[2].y);
			transformInv /= glm::dot(transformInv, transformInv);
			delta.y = glm::dot(delta0, transformInv);
			transformInv = glm::vec3(ellipsoid[0].z, ellipsoid[1].z, ellipsoid[2].z);
			transformInv /= glm::dot(transformInv, transformInv);
			delta.z = glm::dot(delta0, transformInv);
			X[i] +=  delta ;
			V[i] = glm::vec3(0);
		} 
	}
}

void OnIdle() {

	//Fixed time stepping + rendering at different fps
	if ( accumulator >= timeStep )
    {
        StepPhysics(timeStep); 
        accumulator -= timeStep;
    }
	glutPostRedisplay(); 
}

void StepPhysics(float dt) {

	ComputeForces();
	IntegrateImplicit(dt); 
	EllipsoidCollision(); //此处的collision也并未消耗很多时间
	ApplyProvotDynamicInverse();
}

void main(int argc, char** argv) {


	

	//OMP_PARALLEL
	//{
	//#pragma omp for
	//	cout << "Hello" << ", I am Thread " << omp_get_thread_num() << endl; //支持多个for并行！
	//}
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT Cloth Demo [Implicit Integration-Baraff & Witkin's Model]"); //窗口标题

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);
	glutCloseFunc(OnShutdown);

	glewInit();
	InitGL();
	glutMainLoop();
}
