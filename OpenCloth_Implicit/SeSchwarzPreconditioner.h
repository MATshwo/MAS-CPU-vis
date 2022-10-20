/////////////////////////////////////////////////////////////////////////////////////////////
//  Copyright (C) 2002 - 2022,
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions
//  are met:
//     1. Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//     2. Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//     3. The names of its contributors may not be used to endorse or promote
//        products derived from this software without specific prior written
//        permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
//  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
//  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include "SeCsr.h" //csr
#include "SeMorton.h" //莫顿码
#include "SeAabbSimd.h" //aabbximd
#include "SeCollisionElements.h" //碰撞


SE_NAMESPACE_BEGIN

class SeSchwarzPreconditioner
{

public:
	
	//==== input data
	//初始指向一个空指针
	const SeVec3fSimd* m_positions;		// the mesh nodal positions 网格节点位置 sevec3特殊类型，在vectorsimd头文件中定义

	//==== input mesh topology data for collision computation 输入网格拓扑数据用于碰撞计算

	const Int4* m_edges;				// indices of the two adjacent vertices and the two opposite vertices of edges 输入两个相邻顶点和相对顶点的边的索引，四维向量，每一维都是整数
	const Int4* m_faces;				// indices of the three adjacent vertices of faces, the fourth is useless 三角形三个相邻顶点，第四个位置无效
	
	const SeCsr<int>* m_neighbours;	// the adjacent information of vertices stored in a csr format 用csr格式存储顶点的邻接矩阵

public:
	SeSchwarzPreconditioner() {};
	//==== call before time integration once a frame 时间积分前调用一次  每一帧的时间积分
	void AllocatePrecoditioner(int numVerts, int numEdges, int numFaces);

	//==== call before PCG iteration loop PCG迭代循环前调用
	void PreparePreconditioner(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges);

	//void PreparePreconditioner(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges,
	//	const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);


	//==== call during PCG iterations  PCG迭代期间调用
	void Preconditioning(SeVec3fSimd* z, const SeVec3fSimd* residual, int dim); //r初始取b
	std::vector<int>                    m_MapperSortedGetOriginal;
	int m_numVerts = 0;

private:
	// 变量定义 + 函数声明; 要搞明白变量和函数的类型、作用
	 //顶点个数、边个数&三角形个数
	int m_numEdges = 0;
	int m_numFaces = 0;

	int m_frameIndex = 0; //帧索引

	int m_totalSz = 0; ////每一级别簇的个数*size(32)*1.5 ，然后累加的结果
	int m_numLevel = 0; //多级level的级数

	int m_totalNumberClusters; //total簇

	
	SeAabb<SeVec3fSimd>					m_aabb;


	SeArray2D<SeMatrix3f>				m_hessianMapped;  //二阶hessian映射
	std::vector<int>					m_mappedNeighborsNum; //莫顿空间下每个顶点的邻居个数(包括自己)
	std::vector<int>					m_mappedNeighborsNumRemain; //计算过程中使用的remain变量
	SeArray2D<int>						m_mappedNeighbors;  // 重新排序后顶点的邻居索引列表
	SeArray2D<int>						m_mappedNeighborsRemain;

	SeArray2D<int>						m_CoarseSpaceTables; //系数空间表格
	std::vector<int>                    m_prefixOrignal; //初始前缀
	std::vector<unsigned int>           m_fineConnectMask; 
	std::vector<unsigned int>           m_nextConnectMsk;
	std::vector<unsigned int>           m_nextPrefix; //

	std::vector<SeVec3fSimd>			m_mappedPos; //莫顿码映射后的坐标

	std::vector<Int4>					m_coarseTables;
	//goingnext:整数列表,表示当前索引到下一级的索引的映射 m_goingNext[0]=k表示当前level索引为0的点在下一个level的索引为k;长度很长，包含了所有level顶点的映射：
	//注意下一级level编号不再从0开始，需要加上一级所有顶点的个数作为偏移量
	std::vector<int>					m_goingNext; 
	std::vector<int>					m_denseLevel;
	std::vector<Int2>					m_levelSize;  // .x = current level size 当前level的实际顶点个数(不考虑虚设的顶点-不够一个簇强行凑齐的情况)   .y = (prefixed) current level begin index (前缀)当前级别开始索引：当前level之前所有level的顶点的个数不包含当前(而且是包含虚拟顶点的)

	std::vector<SeVec3fSimd>			m_mappedR;
	std::vector<SeVec3fSimd>			m_mappedZ;
	//std::vector<int>                    m_MapperSortedGetOriginal;          // sorted by morton 
	std::vector<int>                    m_mapperOriginalGetSorted;   
	std::vector<SeMorton64>             m_mortonCode;  //莫顿编码

	SeArray2D<SeMatrix3f>				m_hessian32;   //存储最后的hessian系统矩阵 = 不考虑碰撞(初始输入)+两次碰撞附加hessian
	std::vector<SeMatrix3f>             m_additionalHessian32; //由3*3矩阵构成的向量
	std::vector<float>                  m_invSymR; //系统逆矩阵

	
	int m_stencilNum;
	int m_maxStencilNum;
	std::vector<Stencil>				m_stencils;  
	std::vector<Int5>					m_stencilIndexMapped;

private:

	void ComputeLevelNums(int bankSize); //计算级别个数

	void DoAlllocation(); //定位？
	 
	void ComputeTotalAABB(); //求AABB总体

	void ComputeAABB(); //求AABB

	void SpaceSort(); //空间排序(莫顿码排序)

	void FillSortingData();//有序填充数据 

	void DoingSort(); //排序

	void ComputeInverseMapper(); //求逆

	void MapHessianTable();//二阶矩阵映射

	void PrepareCollisionStencils(const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);

	void MapCollisionStencilIndices(); //碰撞的一些处理

	void ReorderRealtime(); //实时重新排序？



	void BuildConnectMaskL0(); //建立映射L0

	void BuildCollisionConnection(unsigned int* pConnect, const int* pCoarseTable);

	void PreparePrefixSumL0();//

	void BuildLevel1(); //L1

	void BuildConnectMaskLx(int level); //建立LX

	void NextLevelCluster(int level); //下一级的簇

	void PrefixSumLx(int level); //sum

	void ComputeNextLevel(int level);//计算下一级

	void TotalNodes();//总体节点数

	void AggregationKernel();//加速内核

	void AdditionalSchwarzHessian2(SeMatrix3f hessian, std::vector<SeMatrix3f>& pAdditionalHessian, SeArray2D<SeMatrix3f>& pDenseHessian, int v1, int v2, const std::vector<int>& pGoingNext, int nLevel, int vNum);

	void PrepareCollisionHessian(); //碰撞准备

	void PrepareHessian(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges); //二阶矩阵准备
	
	void LDLtInverse512();//LDLT 

	void BuildResidualHierarchy(const SeVec3fSimd* m_cgResidual);

	void SchwarzLocalXSym();

	void CollectFinalZ(SeVec3fSimd* m_cgZ); //最后收集
};

SE_NAMESPACE_END