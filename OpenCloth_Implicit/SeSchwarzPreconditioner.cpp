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

#include <algorithm>
#include <unordered_map>
#include<iostream>
//本代码对应SeSchwarzPreconditioner头文件中的类&函数定义

#include "SeSchwarzPreconditioner.h" //施瓦茨预处理器
#include "SeIntrinsic.h" //Se intrinsic内在的


SE_USING_NAMESPACE

void SeSchwarzPreconditioner::AllocatePrecoditioner(int numVerts, int numEdges, int numFaces)
{
	//预处理分配：分配内存之类的
	m_numVerts = numVerts;
	m_numEdges = numEdges;
	m_numFaces = numFaces;

	if (m_frameIndex == 0)
	{
		DoAlllocation();
	}

	if (m_frameIndex % 17 != 0)
	{
		return;
	}

	//==== Reorder 

	ComputeTotalAABB();

	SpaceSort();

	ComputeInverseMapper();

	MapHessianTable();

	m_frameIndex++;
}

void SeSchwarzPreconditioner::PreparePreconditioner
(
	const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges
	//const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges,
	//const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets,  //碰撞元素中定义的三个结构体类型
	//unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts
)
{
	m_mappedNeighborsRemain = m_mappedNeighbors;
	m_mappedNeighborsNumRemain = m_mappedNeighborsNum;
	//PrepareCollisionStencils(efSets, eeSets, vfSets, efCounts, eeCounts, vfCounts); 
	ReorderRealtime();
	const int blockNum = m_totalSz / 32; 
	Utility::MemsetZero(m_hessian32.Ptr(), m_hessian32.Size()); 
	Utility::MemsetZero(m_additionalHessian32); 

	//PrepareCollisionHessian(); //碰撞附加
	PrepareHessian(diagonal, csrOffDiagonals, csrRanges); //对角+压缩对角+压缩范围
	LDLtInverse512(); //LDLt求逆：这个是对矩阵求逆的运算
}

void SeSchwarzPreconditioner::Preconditioning(SeVec3fSimd* z, const SeVec3fSimd* residual, int dim)
{
	//预处理——分级
	Utility::MemsetZero(m_mappedZ); 
	Utility::MemsetZero(m_mappedR);

	BuildResidualHierarchy(residual); 

	SchwarzLocalXSym(); 
	 
	CollectFinalZ(z); 
}

void SeSchwarzPreconditioner::ComputeLevelNums(int bankSize)
{
	constexpr float SizeRatio = 1.5f; //banksize=32，簇规模

	//====参数声明，并初始化；

	int totalSz = 0;
	int nLevel = 1; //level1出发

	//e.g顶点数=64，计算得到的结果就是2*32；如果=65，结果就是3*32；按照bank_size的规模调整格式
	int levelSz = (m_numVerts + bankSize - 1) / bankSize * bankSize; //这里先做整除得到倍数再乘以bank_size得到32整数倍的规模，因为原始顶点个数可能不是32的倍数，导致无法按照簇进行划分
	totalSz += levelSz; //总共的规模=levelSz

	//首层计算后，得到levelSz = k*32,k表示当前层的簇团个数,32对应簇的规模;
	while (levelSz > 32) // 32表示超级节点规模(一个超级节点包含的节点数,与banksize=32非同一个意思)
	{
		//node num of next level,nodesize=32,好处在于下一级的一个节点对应上级的一个簇
		levelSz /= 32;
		nLevel++;
		levelSz = (levelSz + bankSize - 1) / bankSize * bankSize; //此时新的level_size = 1*32<=32,循环相应结束,3个超级节点只能构成一个簇，所以size=1*32;
		//要注意这里（(levelSz + bankSize - 1) / bankSize）表示当前层节点能构成的簇团数；
		// (levelSz + bankSize - 1) / bankSize * bankSize - 表示当前层的size，因为顶点可能不是32的倍数，这里的size是对原始size进行调整的结果，不涉及下一层level参数的信息，只是对当前层的调整size
		totalSz += levelSz; 
	}

	m_numLevel = nLevel; //total level 
	m_totalSz = totalSz * SizeRatio; //Totalsize = totalsize*1.5?
}

void SeSchwarzPreconditioner::DoAlllocation()
{
	// do all location
	constexpr int bankSize = 32;	//==== size of diagonal block 对角块的规模：这里32指的是顶点个数不是矩阵行的个数

	ComputeLevelNums(bankSize);

	//====一系列初始化resize操作

	m_MapperSortedGetOriginal.resize(m_numVerts); //整数向量:将原始节点排序后的索引
	m_mapperOriginalGetSorted.resize(m_numVerts); //整数向量:莫顿空间的序列到原始节点序列的逆映射
	m_mortonCode.resize(m_numVerts); //莫顿码类型的向量：存储莫顿编码，一个三维顶点对应唯一的莫顿编码

	m_hessian32.Resize(bankSize, m_totalSz); //指定行列：32*(total_k*1.5*32), total_k是原始计算的所有簇(包含level0-3所有的)的个数
	m_mappedZ.resize(m_totalSz); //map z total_size = total_k*32
	m_mappedR.resize(m_totalSz); //map r

	m_mappedPos.resize(m_numVerts); //顶点位置映射,规模同顶点

	m_denseLevel.resize(m_numVerts); //densellevel
	m_goingNext.resize(m_numLevel * m_numVerts); //级数*顶点数
	m_coarseTables.resize(m_numVerts); //粗糙表格记录

	m_prefixOrignal.resize(m_numVerts);//原始前缀
	m_fineConnectMask.resize(m_numVerts); //connect mask 
	m_nextConnectMsk.resize(m_numVerts); 
	m_nextPrefix.resize(m_numVerts);


	int triSz = (1 + 96) * 96 / 2 + 16 * 3; //存储逆运算矩阵4704
	const int blockNum = m_totalSz / 32; //block_num？ 由1.5total_size计算出的blocknum = 原始num*1.5
	m_invSymR.resize(triSz * blockNum); //系统矩阵的逆：
	m_additionalHessian32.resize(m_totalSz + 1);
	//附加hessian阵,大小=totalSz+1,32只是顶点规模,而每个顶点又是由3个坐标构成的
	//m_additionalHessian32是由totalSz+1个3*3对角阵构成的向量

	//====

	m_levelSize.resize(m_numLevel + 1); //总level+1,从第0层,存储最原始的数据
	m_CoarseSpaceTables.Resize(m_numLevel, m_numVerts); //粗糙空间表：不同level和顶点个数：0级对应初始状态

	int maxNeighbours = 0;

	for (int vId = 0; vId < m_numVerts; ++vId)
	{
		maxNeighbours = Math::Max(m_neighbours->Size(vId) + 1, maxNeighbours);//遍历所有节点，得到单个顶点最多相邻的顶点个数+1(包含自身？)
	}

	//一系列resize操作
	m_mappedNeighborsNum.resize(m_numVerts);
	m_mappedNeighborsNumRemain.resize(m_numVerts);//remain

	//向量用小写resize,矩阵类型用Resize
	m_mappedNeighbors.Resize(maxNeighbours, m_numVerts); //存储重排序后的邻居索引:[i,j]存储第j个顶点的第i个邻居下标
	m_mappedNeighborsRemain.Resize(maxNeighbours, m_numVerts);

	const int maxCollisionPerVert = 32; //最大碰撞个数：只考虑域内碰撞
	m_maxStencilNum = m_numVerts * maxCollisionPerVert;
	m_stencils.resize(m_maxStencilNum);
	m_stencilNum = 0;
}

void SeSchwarzPreconditioner::ComputeTotalAABB()
{
	m_aabb = SeAabb<SeVec3fSimd>(); //初始化一个包围盒

	ComputeAABB(); // TODO: omp reduce
	//printf("BVH total end!\n");
}


void SeSchwarzPreconditioner::ComputeAABB()//轴对齐包围盒求解 axis alien bound box
{

//MSVC compiler is stuck with OpenMP version 2.0, and unfortunately for you, reduction(max:) was only introduced with version 3.1 of the OpenMP C/C++ standard (that was in September 2011) 

//#pragma omp parallel for reduction(aabbAdd:m_aabb)
	//printf("BVH action!\n");
	for (int vid = 0; vid < m_numVerts; ++vid) // need omp custom reduction operation
	{
		//printf("inital lower : %f\n", m_aabb.Upper.x);
		m_aabb += m_positions[vid];//遍历顶点，得到包含所有顶点的一个包围盒

	}
	//printf("BVH end!\n");
}

void SeSchwarzPreconditioner::SpaceSort()
{
	FillSortingData(); //空间排序：结合莫顿码和顶点坐标排序
	DoingSort();
}

void SeSchwarzPreconditioner::FillSortingData()
{
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			SeVec3fSimd temp = (m_positions[vid] - m_aabb.Lower) / m_aabb.Extent(); //临时变量：将顶点坐标缩放至[0,1]之间

			SeMorton64 mt;

			mt.Encode(temp.x, temp.y, temp.z); //得到当前顶点的莫顿码

			m_mortonCode[vid] = mt; //莫顿码结果存储

			m_MapperSortedGetOriginal[vid] = vid; //原始的排序就是顶点索引的排序，之后还需要调整，初始[0,1,2,3,4,5,6]
		}
}


void SeSchwarzPreconditioner::DoingSort()
{
	//通过莫顿码得到原始索引：第index个莫顿顶点的value是原始索引
	//==== sort m_MapperSortedGetOriginal by m_mortonCode1; 通过莫顿码将原始顶点重新排序(从小到大)
	//[0,1,2,3,4,5,6] -> [1,3,0,6,5,4,2]，表示顶点i排序后的结果，[1,3,0,6,5,4,2]; （value,index）原始第value个元素是莫顿码的第index个元素
	//实际操作：根据排序后的列表[1,3,0,6...]依次找到原始编号为1,3,0,6的顶点按顺序重新编码为0,1,2...
	std::sort(m_MapperSortedGetOriginal.begin(), m_MapperSortedGetOriginal.end(), [&](int i, int j) { return m_mortonCode[i] < m_mortonCode[j]; }); //前两个参数指定范围，第三个指定排序的规则
}

void SeSchwarzPreconditioner::ComputeInverseMapper()
{
	//第index个原始顶点的value是莫顿排序后的索引
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			int originalIndex = m_MapperSortedGetOriginal[vid];

			m_mapperOriginalGetSorted[originalIndex] = vid; //与上边对应的就是[2,0,6,1,5,4,3],表示莫顿码中的第value个元素是原始第index个元素（value,index）逆映射
		}
}


void SeSchwarzPreconditioner::MapHessianTable()
{
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			const auto  originalIndex = m_MapperSortedGetOriginal[vid]; //得到原始按莫顿规则从小到大的index

			m_mappedPos[vid] = m_positions[originalIndex]; //将顶点坐标赋值给莫顿顶点位置

			//====

			int numNeighbors = m_neighbours->Size(originalIndex); //邻居个数赋值

			const int* neighbors = m_neighbours->IdxPtr(originalIndex); //当前顶点所有邻居的列索引列表

			m_mappedNeighborsNum[vid] = numNeighbors + 1;	// include self 该列表存储包含自身的邻居顶点个数

			m_mappedNeighbors[0][vid] = vid;				// include self 存储包含自身的邻居顶点的索引,第vid个顶点的第0个邻居是自己

			for (int k = 1; k < numNeighbors + 1; ++k)
			{
				unsigned int originalNeighbors = neighbors[k - 1];//m_neighbors[k][originalIndex];

				m_mappedNeighbors[k][vid] = m_mapperOriginalGetSorted[originalNeighbors];
			}
		}
}

void SeSchwarzPreconditioner::MapCollisionStencilIndices()
{
	//碰撞组件与莫顿空间的映射
	m_stencilIndexMapped.resize(m_stencilNum);

	OMP_PARALLEL_FOR

		for (int i = 0; i < m_stencilNum; i++)
		{
			const auto& s = m_stencils[i];

			for (int vi = 0; vi < s.verextNumPerStencil; vi++)
			{
				//按莫顿码重新对碰撞体编号:碰撞体1(vf结构)：原始index是[1,3,5,6] -> 到莫顿空间就是 [4,7,2,8]
				//m_stencilIndexMapped[i]表示在莫顿空间中第i个碰撞结果的index结构
				m_stencilIndexMapped[i][vi] = m_mapperOriginalGetSorted[s.index[vi]]; //转换映射 
			}
		}
}

void SeSchwarzPreconditioner::PrepareCollisionStencils(const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts)
{
	//参考文献：PSCC:三角形之间的ee&vf一共有15种情形
	/*实例:4个顶点, 5条边, 2个面的简单案例,{0,1,2,3};efset=[{00}{10}{20}{21}{31}{41}];eeset=[{01}{02}{03}{12}{14}{23}{24}{34}];vfset=[{00}{10}{20}{11}{21}{31}],efcount=[],eecount=[],vfcount=[]*/

	//三类碰撞：线面\点面\线线 
	int efNum = efCounts[m_numEdges]; //count计算,累加列表,num对应个数
	int eeNum = eeCounts[m_numEdges];
	int vfNum = vfCounts[m_numVerts];

	//所有可能碰撞的情形
	int totalStencilNum = efNum + eeNum + vfNum;

	if (totalStencilNum > m_maxStencilNum)
	{
		totalStencilNum = m_maxStencilNum;
		printf("stencil size %d exceed max stencils num  %d\n", totalStencilNum, m_maxStencilNum);
	}

	m_stencilNum = 0;

	OMP_PARALLEL_FOR

		for (int i = 0; i < totalStencilNum; i++)
		{
			Stencil s;

			//创建碰撞组件的列表: 顶点,一个索引;边,两个索引(端点);面,三个索引(顶点);
			//e.g:{ee:每条边需要两个顶点索引确定,所以verextNumPerStencil=2+2}；{ef:,所以verextNumPerStencil=2+3}
			if (i < efNum) //process ef
			{
				const auto& pair = efSets[i]; //第i个e-f连接处

				if (pair.m_eId < 0 || pair.m_fId < 0) { continue; }

				const auto& edge = m_edges[pair.m_eId];
				const auto& face = m_faces[pair.m_fId];

				s.verextNumPerStencil = 5;
				s.vertexNumOfFirstPrimitive = 2;

				s.index[0] = edge[0]; //edge:输入数据，一条边由两个端点构成
				s.index[1] = edge[1];
				s.index[2] = face[0]; //三角面
				s.index[3] = face[1];
				s.index[4] = face[2];

				//权重值与重心的坐标有关
				s.weight[0] = pair.m_bary[0]; //bary重心,三个分量 bary.x,.y,.z
				s.weight[1] = 1.f - pair.m_bary[0];
				s.weight[2] = -pair.m_bary[1];
				s.weight[3] = -pair.m_bary[2];
				s.weight[4] = -(1.f - pair.m_bary[1] - pair.m_bary[2]);

				s.direction = pair.m_normal; //方向与法线一致

				s.stiff = pair.stiff; //刚度系数

			}
			else if (i < efNum + eeNum) // process ee
			{
				const auto& pair = eeSets[i];

				if (pair.m_eId1 < 0 || pair.m_eId0 < 0) { continue; }

				s.verextNumPerStencil = 4;
				s.vertexNumOfFirstPrimitive = 2;

				const auto& edge0 = m_edges[pair.m_eId0];
				const auto& edge1 = m_edges[pair.m_eId1];

				s.index[0] = edge0[0];
				s.index[1] = edge0[1];
				s.index[2] = edge1[0];
				s.index[3] = edge1[1];

				s.weight[0] = pair.m_bary[0];
				s.weight[1] = 1.f - pair.m_bary[0];
				s.weight[2] = -pair.m_bary[1];
				s.weight[3] = -(1.f - pair.m_bary[1]);

				s.direction = pair.m_normal;

				s.stiff = pair.stiff;
			}
			else if (i < totalStencilNum) //process vf 
			{
				const auto& pair = vfSets[i];

				if (pair.m_vId < 0 || pair.m_fId < 0) { continue; }

				s.verextNumPerStencil = 4; 
				s.vertexNumOfFirstPrimitive = 3; //原始顶点num:面-顶点,面为原始 = 3
				
				const auto& face = m_faces[pair.m_fId];
				
				s.index[0] = face[0];
				s.index[1] = face[1];
				s.index[2] = face[2];
				s.index[3] = pair.m_vId;

				s.weight[0] = -pair.m_bary[0];
				s.weight[1] = -pair.m_bary[1];
				s.weight[2] = -(1.f - pair.m_bary[2]);
				s.weight[3] = 1.f;

				s.direction = pair.m_normal;

				s.stiff = pair.stiff;
			}

			int stencilNum = Intrinsic::AtomicAdd(&m_stencilNum, 1); // 等价于执行stencilNum = m_stencilNum;m_stencilNum+=1;

			m_stencils[stencilNum] = s; //将碰撞结构存储到m_stencils中
		}

	MapCollisionStencilIndices(); //1
}

void SeSchwarzPreconditioner::ReorderRealtime()
{
	//实时排序
	Utility::MemsetZero(m_levelSize);

	BuildConnectMaskL0(); //从level-0开始

	BuildCollisionConnection(m_fineConnectMask.data(), nullptr); //2

	PreparePrefixSumL0(); 

	BuildLevel1(); 

	for (int level = 1; level < m_numLevel; level++)
	{
		//遍历level
		Utility::MemsetZero(m_nextConnectMsk);

		BuildConnectMaskLx(level); //建立mask

		BuildCollisionConnection(m_nextConnectMsk.data(), m_CoarseSpaceTables[level - 1]); //碰撞联系

		NextLevelCluster(level); //

		PrefixSumLx(level); //

		ComputeNextLevel(level); //计算下一个level
	}

	TotalNodes(); //节点总数

	AggregationKernel(); //加速核
}

void SeSchwarzPreconditioner::BuildConnectMaskL0()
{
	//首先明确建立的connetmask是何含义? 指代各级顶点之间的连接关系:相同簇内的连接信息和不同簇内的连接信息(随着level深入,还要考虑顶点合并信息更新当前level的顶点连接信息)
	//mask记录了簇内顶点之间的相邻信息，不考虑簇外！

	//对每个原始顶点在level0分簇，并更新其mask和neighbor信息
	//L0的联系mask? 这部分应该对应对原始域？
	const unsigned int warpDim = 32; //扭曲维度 = banksize
	int nWarps = (m_numVerts + warpDim - 1) / warpDim; //簇团数 warp:包,这里就理解为簇团数

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < nWarps; ++warpIdx) 
		{	// warp-id:第几个簇 0-nwarp
			for (int laneIdx = 0; laneIdx < warpDim; ++laneIdx)
			{
				//laneid:簇内第几个顶点0-31
				int vId = warpIdx * warpDim + laneIdx; // vID: 0-nWarps*warpdim/numvert

				if (vId >= m_numVerts)
				{
					//超过顶点总数跳出循环，当顶点个数不是32的倍数有可能发生
					break;
				}
				int numNeighbors = m_mappedNeighborsNumRemain[vId]; //当前顶点的邻居个数,mapped说明在莫顿空间上操作
				/*1
				2
				4
				8
				16
				32
				64
				128
				256
				512
				1024
				2048
				4096
				8192
				16384
				32768
				65536
				131072
				262144
				524288
				1048576
				2097152
				4194304
				8388608
				16777216
				33554432
				67108864
				134217728
				268435456
				536870912
				1073741824
				2147483648 2^0 - 2^31 内部循环一次得到初始的mask形式
				*/
				unsigned int connetMsk = 1U << laneIdx;	// include self 1U-32位无符号整型,左移idx个二进制位

				int nk = 0;

				for (int k = 0; k < numNeighbors; k++)
				{
					int vIdConnected = m_mappedNeighborsRemain[k][vId]; //得到第k个邻居的索引

					int warpIdxConnected = vIdConnected / warpDim;//得到簇的索引，表示vID顶点与簇连接的信息:e.g vID=1,与索引34的顶点相连，则vID与簇1相连接

					if (warpIdx == warpIdxConnected) // in the same warp 若当前顶点与其邻居在同一个簇内
					{
						unsigned int laneIdxConnected = vIdConnected % warpDim; //取余数

						//更新mask,按位或后赋值,按道理应该等于0
						connetMsk |= (1U << laneIdxConnected);  //相同簇内若有邻居:mask = (1U << laneIdx) | (1U << laneIdxConnected),结果:有一个是1的位=1;都是0的位=0,数值增大
						//再举个例子:以编号为0的质点为例，初始mask=1U<<0 = 0001 = 1;簇大小是4(二进制位只需要4位即可)，相邻点索引{0,2,6},{0,2}与0在相同簇,进入循环,0%4=0,2%4=2,最后的mask = 0001|0100 = 0101 
						//观察顶点的connetmask可以直接观察到簇内(只能是簇内)顶点的相邻信息
					}
					else
					{
						m_mappedNeighborsRemain[nk][vId] = vIdConnected; //当前顶点存在连接点跨越不同簇，则作为remain邻居(留下的邻居),并将索引重新储存到remainneighbor中
						nk++;
					}
				}
				m_mappedNeighborsNumRemain[vId] = nk; //更新当前节点的邻居信息，比之前显著减少,减少个数 = 与当前顶点在同一簇内且相连接的顶点

				// distance culling 
				//const SeVec3fSimd& pos = m_mappedPos[vId];
				//for (int it = 0; it < Intrinsic::Popcount32(activeMsk); it++)
				//{
				//	const float th = AABB_TH * AABB_TH;
				//	CuFloat3 po(0.0f);
				//	po.x = Intrinsic::Shuffle(activeMsk, pos.x, it);
				//	po.y = Intrinsic::Shuffle(activeMsk, pos.y, it);
				//	po.z = Intrinsic::Shuffle(activeMsk, pos.z, it);
				//	if (Math::SqrLength(po - pos) < th)
				//	{
				//		connetMsk |= (1U << it);
				//	}
				//}
				
				//finemask记录了每个顶点与所在簇其他顶点的连接信息包含自身
				//m_mappedNeighborsRemain[:][vertexindex]每个顶点与其他簇中顶点的连接信息:remainneighbor不包含自身
				//将簇内和簇外的连接视作两种不同的情况加以区分 - 二者连接数量相加=原始邻居个数
				m_fineConnectMask[vId] = connetMsk; //存储每个节点的mask信息,fine存储最后全部的mask结果,长度同Vertnum
			}
		}
}

void SeSchwarzPreconditioner::BuildCollisionConnection(unsigned int* pConnect, const int* pCoarseTable)
{
	//建立碰撞联系:碰撞部分暂时放下
	const int bank = 32;

	OMP_PARALLEL_FOR

		for (int i = 0; i < m_stencilNum; i++)
		{
			const auto& s = m_stencils[i];
			Int5 idx = m_stencilIndexMapped[i];
			unsigned int connMsk[5] = {};
			if (pCoarseTable)
			{
				for (int it = 0; it < s.verextNumPerStencil; it++)
				{
					idx[it] = pCoarseTable[idx[it]];
				}
			}
			for (int ita = 0; ita < s.verextNumPerStencil; ita++)
			{
				for (int itb = ita + 1; itb < s.verextNumPerStencil; itb++)
				{
					unsigned int myID = idx[ita];
					unsigned int otID = idx[itb];
					if (myID == otID) // redundant ? 多余的？
					{
						continue;
					}
					if ((myID / bank == otID / bank))
					{
						if (ita < s.vertexNumOfFirstPrimitive && itb >= s.vertexNumOfFirstPrimitive)
						{
							connMsk[ita] |= (1U << (otID % bank));
							connMsk[itb] |= (1U << (myID % bank));
						}
					}
				}
			}
			for (int it = 0; it < s.verextNumPerStencil; it++)
			{
				if (connMsk[it])
				{
					Intrinsic::AtomicOr(&pConnect[idx[it]], connMsk[it]);
				}
			}
		}
}

void SeSchwarzPreconditioner::PreparePrefixSumL0()
{
	//此处实现了连接顶点之间的跳转:分别当前节点跳转到与其连接的每一个顶点，并更新mask,保证连接的顶点构成一个闭环
	
	//常规: 簇规模(dim),簇团数
	int warpDim = 32;
	int nWarps = (m_numVerts + warpDim - 1) / warpDim;

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < nWarps; ++warpIdx)
		{
			unsigned int cacheMask[32]; //临时缓存变量:存储mask信息;长度32,外部循环一次更新一次

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vId = laneIdx + warpIdx * warpDim; //当前顶点索引

				if (vId >= m_numVerts) { break; }

				cacheMask[laneIdx] = m_fineConnectMask[vId]; //加载数据到cache or cache列表赋值
			}

			Intrinsic::SyncWarp(); //CUDA线程同步

			int prefixSum = 0; //计数变量

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * warpDim; //真实顶点的索引

				if (vid >= m_numVerts) { break; }

				unsigned int connetMsk = cacheMask[laneIdx];//提取cache中数据——connet

				unsigned int visited = (1U << laneIdx); //按位左移: 将1U左移landIdx位,计算的结果=1U*(2^(laneIdx))

				//unsigned int = -1,表示无符号整型的最大值 = 11111111111111111111111111111111,此时说明当前顶点与簇内所有顶点连接
				
				//该循环为了保证不存在单向连接，所有连接的顶点之间保证连通——对mask再次更新
				while (connetMsk != -1) 
				{
					//遍历connetmsk中非0的每一位,visited记录访问状态,todo指示变量，表示是否遍历完毕！
					unsigned int todo = visited ^ connetMsk; // 对二进制的每一位执行相同位=0;不同位=1

					if (!todo)
					{
						//todo = 0跳出循环:即visited == connetMsk 每一位都相等
						break;
					}

					unsigned int nextVisit = Intrinsic::Ffs(todo) - 1; //返回todo首个非零位的索引

					visited |= (1U << nextVisit); //跳转到下一个连接的点，更新visit状态

					connetMsk |= cacheMask[nextVisit]; //更新mask,为了保证所有域内连接点是连通的
				}
				m_fineConnectMask[vid] = connetMsk; //存储更新后的mask信息

				//__popcount(temp)计算temp二进制位中=1的个数; &遇0则0;|遇1则1(长得也像1)
				//LanemaskLt:(1U << laneId) - 1 = 2^(laneId)-1;作用:计数connetMsk前laneIdx二进制位中=1的个数
				//此处执行的功能存疑？计算返回的结果恰好是二进制位非零位索引+1,有可以直接调用的函数为什么还要用这个？
				unsigned int electedPrefix = Intrinsic::Popcount32(connetMsk & Intrinsic::LanemaskLt(laneIdx)); 
				if (electedPrefix == 0)
				{
					//=0表示与前id个顶点都不连接,sum+1
					prefixSum++;
				}
			}
			m_prefixOrignal[warpIdx] = prefixSum; //首个连接的顶点的位置？
		}
}

void SeSchwarzPreconditioner::BuildLevel1()
{
	constexpr unsigned int blockDim = 1024;//32*32

	int nBlocks = (m_numVerts + blockDim - 1) / blockDim;//level1顶点数 = level0/32，level1的簇团数= level0/32/32

	OMP_PARALLEL_FOR
	// 结构:for{ for1{ for11};for2{};for3{for31 for32 for33 for34}   }

		//遍历level1的簇:block>(super node)>(node)
		for (int blockIdx = 0; blockIdx < nBlocks; ++blockIdx) //遍历簇
		{
			unsigned int warpPrefix[32];//临时变量，用于保存每个簇内的prefix信息,长度32,每次大循环更新

			unsigned int theBlockGlobalPrefix = 0; //block全局前缀

			unsigned int globalWarpOffset = blockIdx * 32; //簇偏移:一个簇偏移32个长度

			//还不明确具体含义,得到theBlockGlobalPrefix
			for (int warpIdx = 0; warpIdx < blockIdx; ++warpIdx) //
			{
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;
					theBlockGlobalPrefix += m_prefixOrignal[threadIdx]; //叠加,而且有重新叠加的行为，这个数量非常大
				}
			}

			//将prefix信息拷贝到临时列表中
			for (int warpIdx = 0; warpIdx < 32; ++warpIdx)
			{
				warpPrefix[warpIdx] = m_prefixOrignal[globalWarpOffset + warpIdx];
			}

			Intrinsic::SyncBlock();

			//计算coarseTable & goingMap
			for (unsigned int warpIdx = 0; warpIdx < 32; ++warpIdx) //32-簇内遍历
			{
				unsigned int theWarpGlobalPrefix = theBlockGlobalPrefix;

				for (unsigned int prevWarpIdx = 0; prevWarpIdx < warpIdx; ++prevWarpIdx)
				{
					theWarpGlobalPrefix += warpPrefix[prevWarpIdx]; //再次叠加
				}

				//计算electedmask
				unsigned int electedMask = 0;
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32; //level1下超级顶点的编号(thread-supernode)

					int vid = threadIdx + blockIdx * blockDim; //initial vertex index

					if (vid >= m_numVerts) { break; }


					if (vid == m_numVerts - 1) //最后一个顶点
					{
						m_levelSize[1].x = theWarpGlobalPrefix + warpPrefix[warpIdx];
						m_levelSize[1].y = (m_numVerts + 31) / 32 * 32;
					}

					unsigned int connMsk = m_fineConnectMask[vid];

					unsigned int electedPrefix = Intrinsic::Popcount32(connMsk & Intrinsic::LanemaskLt(laneIdx));

					if (electedPrefix == 0) 
					{
						electedMask |= 1U << laneIdx; 
					}
					//unsigned int electedMask = Intrinsic::BallotSync(-1, electedPrefix == 0);
				}

				//计算lanemake[],长度32
				unsigned int lanePrefix[32];
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= m_numVerts) { break; }


					lanePrefix[laneIdx] = Intrinsic::Popcount32(electedMask & Intrinsic::LanemaskLt(laneIdx));

					lanePrefix[laneIdx] += theWarpGlobalPrefix;
				}

				//计算coarseTable & goingnext
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= m_numVerts) { break; }


					unsigned int connMsk = m_fineConnectMask[vid];

					unsigned int elected_lane = Intrinsic::Ffs(connMsk) - 1;

					unsigned int theLanePrefix = lanePrefix[elected_lane];

					m_CoarseSpaceTables[0][vid] = theLanePrefix;

					m_goingNext[vid] = theLanePrefix + (m_numVerts + 31) / 32 * 32;  //goingNext构造了每一级的质点到下一级的映射(注意是包含所有level的所有顶点！)
				}
			}
		}
}


void SeSchwarzPreconditioner::BuildConnectMaskLx(int level)
{
	//建立不同层之间的联系
	constexpr unsigned int warpDim = 32;

	int nWarps = (m_numVerts + warpDim - 1) / warpDim;

	const unsigned int bank = 32;

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < nWarps; ++warpIdx)
		{
			unsigned int prefixMsk[32];
			unsigned int connectMsk[32];

			for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= m_numVerts) { break; }

				prefixMsk[laneIdx] = m_fineConnectMask[vid];

				unsigned int coarseVid = m_CoarseSpaceTables[level - 1][vid];

				connectMsk[laneIdx] = 0;// 1U << (coarseVid % bank);	//==== include self

				unsigned int kn = m_mappedNeighborsNumRemain[vid];

				unsigned int nk = 0;

				for (unsigned int k = 0; k < kn; k++)
				{
					unsigned int connect = m_mappedNeighborsRemain[k][vid];

					unsigned int coarseConnect = m_CoarseSpaceTables[level - 1][connect];

					if (coarseVid / bank == coarseConnect / bank)
					{
						unsigned int off = coarseConnect % bank;

						connectMsk[laneIdx] |= (1U << off);
					}
					else
					{
						m_mappedNeighborsRemain[nk][vid] = connect;
						nk++;
					}
				}

				m_mappedNeighborsNumRemain[vid] = nk;
			}

			bool isFullWarp = (prefixMsk[0] == -1);

			if (isFullWarp)
			{
				unsigned int connectMskFull = 0;

				for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					connectMskFull |= connectMsk[laneIdx];
				}
				for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					connectMsk[laneIdx] = connectMskFull;
				}
			}
			else
			{
				unsigned int cacheMsk[32];

				for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					cacheMsk[laneIdx] = 0;
				}

				unsigned int electedLane[32];

				for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					electedLane[laneIdx] = Intrinsic::Ffs(prefixMsk[laneIdx]) - 1;

					if (connectMsk[laneIdx])
					{
						//cacheMsk[electedLane[laneIdx]] |= connectMsk[laneIdx];
						Intrinsic::AtomicOr(&cacheMsk[electedLane[laneIdx]], connectMsk[laneIdx]);
					}
				}

				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					connectMsk[laneIdx] = cacheMsk[electedLane[laneIdx]];
				}
			}

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= m_numVerts) { break; }

				unsigned int coarseVid = m_CoarseSpaceTables[level - 1][vid];

				unsigned int electedPrefix = Intrinsic::Popcount32(prefixMsk[laneIdx] & Intrinsic::LanemaskLt(laneIdx));

				if (connectMsk[laneIdx] && electedPrefix == 0)
				{
					Intrinsic::AtomicOr(m_nextConnectMsk.data() + coarseVid, connectMsk[laneIdx]);
				}
			}
		}
}

void SeSchwarzPreconditioner::NextLevelCluster(int level)
{
	//下一级的簇计算
	const int levelNum = m_levelSize[level].x;

	const int warpDim = 32; //

	int numWarps = (m_numVerts + warpDim - 1) / warpDim;

	int maxWarpIdx = (levelNum + warpDim - 1) / warpDim;

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < numWarps; ++warpIdx)
		{
			unsigned int cachedMsk[32];

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= (levelNum + 31) / 32 * 32)
				{
					break;// early culling
				}

				unsigned int connectMsk = (1U << laneIdx);

				if (vid < levelNum)
				{
					connectMsk |= m_nextConnectMsk[vid];// connect to current 
				}

				cachedMsk[laneIdx] = connectMsk;

				if (vid >= levelNum)
				{
					break;
				}
			}

			Intrinsic::SyncWarp();

			unsigned int prefixSum = 0;

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= levelNum || vid >= (levelNum + 31) / 32 * 32)
				{
					break;// early culling
				}

				unsigned int connectMsk = cachedMsk[laneIdx];

				unsigned int visited = (1U << laneIdx);

				while (true)
				{
					unsigned int todo = visited ^ connectMsk;

					if (!todo)
					{
						break;
					}

					unsigned int nextVisist = Intrinsic::Ffs(todo) - 1;

					visited |= (1U << nextVisist);

					connectMsk |= cachedMsk[nextVisist];
				}

				m_nextConnectMsk[vid] = connectMsk;

				unsigned int electedPrefix = Intrinsic::Popcount32(connectMsk & Intrinsic::LanemaskLt(laneIdx));

				if (electedPrefix == 0)
				{
					prefixSum++;
				}
			}

			if (warpIdx < maxWarpIdx)
			{
				m_nextPrefix[warpIdx] = prefixSum;
			}
		}
}

void SeSchwarzPreconditioner::PrefixSumLx(int level)
{
	// 前缀求和
	const int levelNum = m_levelSize[level].x;

	const int levelBegin = m_levelSize[level].y;

	unsigned int blockDim = 1024;//32*32

	int nBlocks = (m_numVerts + blockDim - 1) / blockDim; //这个是下一层的簇的个数

	OMP_PARALLEL_FOR

		for (int blockIdx = 0; blockIdx < nBlocks; ++blockIdx)
		{
			unsigned int warpPrefix[32];

			unsigned int theBlockPrefix = 0;

			unsigned int globalWarpOffset = blockIdx * 32;

			for (int warpIdx = 0; warpIdx < blockIdx; ++warpIdx)
			{
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					unsigned int vid = threadIdx + blockIdx * blockDim;

					if (vid >= (levelNum + blockDim - 1) / blockDim * blockDim)
					{
						break;
					}

					theBlockPrefix += m_nextPrefix[threadIdx];
				}
			}

			for (int warpIdx = 0; warpIdx < 32; ++warpIdx)
			{
				warpPrefix[warpIdx] = m_nextPrefix[globalWarpOffset + warpIdx];
			}

			Intrinsic::SyncBlock();

			for (unsigned int warpIdx = 0; warpIdx < 32; ++warpIdx)
			{
				unsigned int theWarpPrefix = theBlockPrefix;

				for (unsigned int prevWarpIdx = 0; prevWarpIdx < warpIdx; ++prevWarpIdx)
				{
					theWarpPrefix += warpPrefix[prevWarpIdx];
				}

				unsigned int electedMsk = 0;

				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= levelNum)
					{
						break;
					}
					if (vid == levelNum - 1)
					{
						m_levelSize[level + 1].x = theWarpPrefix + warpPrefix[warpIdx];
						m_levelSize[level + 1].y = levelBegin + (levelNum + 31) / 32 * 32;
					}

					unsigned int connMsk = m_nextConnectMsk[vid];

					unsigned int electedPrefix = Intrinsic::Popcount32(connMsk & Intrinsic::LanemaskLt(laneIdx));

					if (electedPrefix == 0)
					{
						electedMsk |= (1U << laneIdx);
					}

					//unsigned int electedMask = Intrinsic::BallotSync(-1, electedPrefix == 0);
				}

				unsigned int lanePrefix[32];

				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					lanePrefix[laneIdx] = Intrinsic::Popcount32(electedMsk & Intrinsic::LanemaskLt(laneIdx));
					lanePrefix[laneIdx] += theWarpPrefix;
				}

				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= levelNum) { break; }

					int electedLane = Intrinsic::Ffs(m_nextConnectMsk[vid]) - 1;

					unsigned int theLanePrefix = lanePrefix[electedLane];

					m_nextConnectMsk[vid] = theLanePrefix;

					m_goingNext[vid + levelBegin] = theLanePrefix + levelBegin + (levelNum + 31) / 32 * 32;
				}
			}
		}
}

void SeSchwarzPreconditioner::ComputeNextLevel(int level)
{
	OMP_PARALLEL_FOR
		//计算下个level
		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			int next = m_CoarseSpaceTables[level - 1][vid];

			m_CoarseSpaceTables[level][vid] = m_nextConnectMsk[next];
		}
}

void SeSchwarzPreconditioner::TotalNodes()
{
	//第i层的总的簇个数
	m_totalNumberClusters = m_levelSize[m_numLevel].y;
	//printf("%d\n", m_totalNumberClusters);
}

void SeSchwarzPreconditioner::AggregationKernel()
{
	const int warpDim = 32;

	int numWarps = (m_numVerts + warpDim - 1) / warpDim; //第一级的cluster数量

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < numWarps; ++warpIdx)
		{
			unsigned int firstInWarp = warpIdx * 32;

			int curr[32];
			int next[32];
			int aggLevel[32];

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= m_numVerts) { break; }

				curr[laneIdx] = vid; //长度32,隶属于簇内的变量

				aggLevel[laneIdx] = m_numLevel - 1;
			}

			Int4 coarseTable[32];

			for (int l = 0; l < m_numLevel - 1; ++l)
			{
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					next[laneIdx] = m_goingNext[curr[laneIdx]];
					//e.g:现在是顶点0,curr[0]=32表示当前簇内第0个顶点是初始第32个顶点;m_goingNext[curr[0]]=66则next[0]=1表示当前簇内第0个顶点是下一个level索引为7的顶点
					//相当于:以laneidx为索引,curr存储了在当前level对应的顶点索引;next存储了映射到下一级的顶点索引
				}

				Intrinsic::SyncWarp();

				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					int vid = laneIdx + warpIdx * 32;

					if (vid >= m_numVerts) { break; }

					if (next[laneIdx] == next[0])
					{
						aggLevel[laneIdx] = Math::Min(l, aggLevel[laneIdx]);
					}

					curr[laneIdx] = next[laneIdx];
					coarseTable[laneIdx][l] = next[laneIdx]; //这里是L不是1
				}
			}


			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= m_numVerts) { break; }

				m_denseLevel[vid] = aggLevel[laneIdx]; //agglevel&denselevel在之后计算中貌似都没有用到;暂不考虑

				m_coarseTables[vid] = coarseTable[laneIdx]; 
				//coarse长度为32;m_coarseTables长度为numvert:记录了每个顶点所有的映射信息
				//猜想：形式类似m_coarseTable = [[65,66,67],[65,66,67],[65,66,67]]...其中索引为1的列表[0,0,0],表示初始第一个顶点分别映射到level1,2,3的索引是65,66,67;
			}
		}
}

void SeSchwarzPreconditioner::AdditionalSchwarzHessian2(SeMatrix3f hessian, std::vector<SeMatrix3f>& pAdditionalHessian, SeArray2D<SeMatrix3f>& pDenseHessian, int v1, int v2, const std::vector<int>& pGoingNext, int nLevel, int vNum)
{
	//本函数的作用：考虑同一域内节点的碰撞
	//1) additional 这个附加是在二者能够合并域所在level映射后的超级顶点，再映射到下一个level的超级顶点(对比下一行又映射了一次)，然后进行赋值！所以这部分不包含最后一层赋值，因为最后一层无法再映射到下一层！
	//2) dense 原始碰撞的两个节点一直循环到他们能够合并的域，则在二者能够合并为一个域所在level的映射后的超级节点列索引处对hessian32进行赋值(记录了原始节点映射到某个level同一个域内超级节点的碰撞)；这是一个附加



	int level = 0;
	unsigned int myID = v1;
	unsigned int otID = v2;
	const unsigned int bank = 32; //簇的大小

	while (myID / bank != otID / bank && level < nLevel)
	{
		//两个点不能划到同一个域跳则到下一级节点直至二者可以合并到同一个域或者到达最大层结束
		myID = pGoingNext[myID];
		otID = pGoingNext[otID];
		level++;
	}

	if (level >= nLevel) //满足第二个循环终止条件则return
		return;
	
	//必有level < nLevel && myID / bank == otID / bank
	//把m_hessian32的3*3对角块先填充好: 每个level都可能附加0-3
	Intrinsic::AtomicAdd(&pDenseHessian[otID % bank][myID], hessian); //hessian32最终的hessian矩阵：在二者当前能够合并的域赋值
	Intrinsic::AtomicAdd(&pDenseHessian[myID % bank][otID], hessian);

	if (level < nLevel - 1)
	{
		//在level1-3附加
		myID = pGoingNext[myID];
		otID = pGoingNext[otID];

		if (myID == otID) //两个顶点在下一个level能够合并为一个节点，比如64，65合并为100
		{
			Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian * 2.0f); //在碰撞顶点能够合并的超级顶点赋值2*hessian
		}
		else
		{
			//参数vnum也没有用
			//此处二者怎么可能不相等呢？想不通 = =；不相等也可以：banksize != supernode size
			//或者说这里的else好像不太对，应该对应外层的那个if?
			std::cout << "Label this running!" << std::endl;
			Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian);
			Intrinsic::AtomicAdd(&pAdditionalHessian[otID], hessian);
		}
	}
	//这里的问题在于：level=3<nlevel=4时，结合上边循环必有 myID / bank == otID / bank,说明他俩在最后一层刚好在一个簇内，但最后一层不再合并顶点；那按照这个代码就会确实最后一层超级顶点的赋值？
	//else 
	//{
	// //这里是我补充的，先注释：else也就是level = nLevel - 1,到了最后一层没有goingnext，且两个节点在一个簇内，但不是同一个节点，这时候需要进行额外赋值
	//	Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian);
	//	Intrinsic::AtomicAdd(&pAdditionalHessian[otID], hessian);
	//}
	//20220915:思考了下原始代码是正确的：考虑此处为什么要叠加，叠加的条件是两个顶点在当前级位于同一个簇，则在两个顶点对应的下一级超级顶点的索引处分别附加一个hessian矩阵；当两个顶点已经位于最后一层时，不会再有下一级，所以也不存在下一级的附加
	//而在当前级别的附加已经赋值到了pdensehessian矩阵中，所以问题不大


}

void SeSchwarzPreconditioner::PrepareCollisionHessian()
{
	//执行碰撞响应：在初始系统矩阵上修改,此处的碰撞响应是在检测到碰撞时间点的上一个时间点额外施加一个外力阻止这个碰撞的发生,表现为叠加在原始系统矩阵上
	//添加排斥势能来接触碰撞——wang2020 repulsion
	OMP_PARALLEL_FOR
		//m_additionalHessian32初始化=0
		for (int i = 0; i < m_stencilNum; i++)
		{
			const auto& s = m_stencils[i];
			auto idx = m_stencilIndexMapped[i]; //使用map后的idx

			Float3 d(s.direction[0], s.direction[1], s.direction[2]); //将s的方向赋值给d

			SeMatrix3f hessian = OuterProduct(d, d * s.stiff); //得到3*3矩阵 ： hessian,初步斥力系数


			for (int it = 0; it < s.verextNumPerStencil; it++)
			{
				
				Intrinsic::AtomicAdd(&m_additionalHessian32[idx[it]], hessian * Math::Square(s.weight[it])); //初次更新m_additionalHessian32[idx[it]]，按莫顿空间的索引更新
			}
			
			for (int ita = 0; ita < s.verextNumPerStencil; ita++)
			{
				for (int itb = ita + 1; itb < s.verextNumPerStencil; itb++)
				{
					//二次附加
					AdditionalSchwarzHessian2(s.weight[ita] * s.weight[itb] * hessian, m_additionalHessian32, m_hessian32, idx[ita], idx[itb], m_goingNext, m_numLevel, m_numVerts);
				}
			}
		}
}

void SeSchwarzPreconditioner::PrepareHessian(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges)
{

	const unsigned int bank = 32;

	int nVC = (m_numVerts + 31) / 32 * 32; //从level1出发？

	OMP_PARALLEL_FOR

		for (int vid = nVC; vid < m_totalNumberClusters; ++vid) //m_totalNumberClusters所有节点的个数: 其中可能包含不存在的节点：因为不够32的被迫凑成32
		{
			auto oldDiagonal = m_additionalHessian32[vid]; //3*3
			int myID = vid;
			Intrinsic::AtomicAdd(&m_hessian32[myID % bank][myID], oldDiagonal); //附加到原始对角阵
			while (true)
			{
				myID = m_goingNext[myID];
				if (myID >= m_totalNumberClusters)
				{
					break;
				}
				Intrinsic::AtomicAdd(&m_hessian32[myID % bank][myID], oldDiagonal);
			}
		}

	int nblock = (m_numVerts + 31) / 32;
	OMP_PARALLEL
	{
		//多线程处理
		std::vector<std::unordered_map<int,SeMatrix3f>> diagTable(m_numLevel);
		
#pragma omp for
		for (int bid = 0; bid < nblock; ++bid)
		{
			for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + bid * 32;
				if (vid >= m_numVerts)
					break;
				const int vidOrig = m_MapperSortedGetOriginal[vid];
				int oldNum = m_mappedNeighborsNum[vid]; //包含自身的neighbor数量

				auto oldDiagonal = diagonal[vidOrig] + m_additionalHessian32[vid];
				m_hessian32[vid % bank][vid] += oldDiagonal;   // self diagonal

				for (int k = 1; k < oldNum; k++)  //index从1开始,不包含自己这个邻居
				{
					const unsigned int neighbor = m_mappedNeighbors[k][vid];
					SeMatrix3f mat = csrOffDiagonals[csrRanges[vidOrig] + k - 1]; //csrrange类似start,存储每个质点首个offdiagnoal在列表csrrange[]中的index = 前i-1个顶点的邻居之和(不含自己)
					int level = 0;
					unsigned int levelSz = m_numVerts;
					unsigned int myID = vid;
					unsigned int otID = neighbor;

					while (myID / bank != otID / bank && level < m_numLevel)
					{
						level++;
						myID = m_goingNext[myID];
						otID = m_goingNext[otID];
					}
					if (level >= m_numLevel)
					{
						continue;
					}
					if (level <= 1) // 按照block分并行，level 1的位置也一定属于本线程
						m_hessian32[otID % bank][myID] += mat;
					else
						Intrinsic::AtomicAdd(&m_hessian32[otID % bank][myID], mat);

					if (level == 0)
						oldDiagonal += mat;
					else if (level + 1 < m_numLevel)
					{
						myID = m_goingNext[myID];
						auto& table = diagTable[level + 1];
						if (table.find(myID) == table.end())
							table[myID] = mat;
						else
							table[myID] += mat;
					}
				}
				if (1 < m_numLevel)
				{
					int myID = m_goingNext[vid];
					m_hessian32[myID % bank][myID] += oldDiagonal;   // self diagonal
					if (2 < m_numLevel)
					{
						myID = m_goingNext[myID];
						auto& table = diagTable[2];
						if (table.find(myID) == table.end())
							table[myID] = oldDiagonal;
						else
							table[myID] += oldDiagonal;
					}
				}
			}
		}

		for (int lv = 2; lv < m_numLevel; lv++)
		{
			const auto& table = diagTable[lv];
			for (auto& pair : table)
			{
				int myID = pair.first;
				Intrinsic::AtomicAdd(&m_hessian32[myID % bank][myID],pair.second);
				if (lv + 1 < m_numLevel)
				{
					myID = m_goingNext[myID];
					auto& nexttable = diagTable[lv + 1];
					if (nexttable.find(myID) == nexttable.end())
						nexttable[myID] = pair.second;
					else
						nexttable[myID] += pair.second;
				}
			}
		}
	}
}

void SeSchwarzPreconditioner::LDLtInverse512()
{
	//LDlt(LDL转置)求逆
	const int triSz = (1 + 96) * 96 / 2 + 16 * 3;

	const int activeblockNum = m_totalNumberClusters / 32; //m_totalNumberClusters全部节点的个数(包含虚拟节点)：activeblockNum构建域的个数

	OMP_PARALLEL_FOR
		//printf("%d\n",CPU_THREAD_NUM);
		for (int block = 0; block < activeblockNum; block++)
		{
			float A[96][96] = {};
			for (int x = 0; x < 32; x++)
			{
				for (int y = 0; y < 32; y++)
				{
					SeMatrix3f temp = m_hessian32[y][x + block * 32]; //m_hessian = (无碰撞系数阵+附加碰撞矩阵)+多级域构建的结果 = final 系数矩阵 = 构建得到域的个数(totalclusters)个32*32矩阵; 

					if (x == y && temp(0, 0) == 0.0f)
					{
						temp = SeMatrix3f::Identity(); //单位阵
					}
					for (int ii = 0; ii < 3; ii++)
					{
						for (int jj = 0; jj < 3; jj++)
						{
							A[x * 3 + ii][y * 3 + jj] = temp(ii, jj); //将A矩阵补齐,右侧添加单位矩阵
						}
					}
				}
			}

			//const int pt = -1;//212;
			//if (block == pt)
			//{
			//	std::ofstream ofs("c:\\test\\MRF.txt");
			//	for (int y = 0; y < 96; y++)
			//	{
			//		for (int x = 0; x < 96; x++)
			//		{
			//			ofs << A[y][x] << " ";
			//		}
			//		ofs << std::endl;
			//	}
			//}

			// 向下消元
#ifdef WIN32
			for (int x = 0; x < 96; x++)
			{
				float diag = A[x][x];
				__m256 line[12];
				for (int it = 0; it < 12; it++)
					line[it] = _mm256_loadu_ps(A[x] + it * 8);

				for (int y = x + 1; y < 96; y++)
				{
					if (A[y][x] == 0.0f)
						continue;
					float r = -A[y][x] / diag;
					__m256 ratio = _mm256_set1_ps(r);
					for (int it = 0; it < 12; it++)
					{
						__m256 temp = _mm256_fmadd_ps(ratio, line[it], _mm256_loadu_ps(A[y] + it * 8)); // A[y] += A[x] * ratio 
						_mm256_storeu_ps(A[y] + it * 8, temp);
					}
					A[y][x] = r;
				}
			}


			// inv diagonal 
			float diagonal[96] = {};
			for (int y = 0; y < 96; y++)
			{
				diagonal[95 - y] = A[y][y];
				A[y][y] = 1.0f;
				for (int x = y + 1; x < Math::Min(96, y + 9); x++)
				{
					A[y][x] = 0.0f;
				}
			}
			for (int it = 0; it < 12; it++)
			{
				__m256 temp = _mm256_div_ps(_mm256_set1_ps(1.0f), _mm256_loadu_ps(diagonal + it * 8));
				_mm256_storeu_ps(diagonal + it * 8, temp);
			}

			int off = block * triSz;
			// output diagonal 
			for (int it = 0; it < 12; it++)
			{
				int lc = 96 - it * 8;
				__m256 acc = _mm256_setzero_ps();
				for (int l = 0; l < lc; l++)
				{
					__m256 a = _mm256_loadu_ps(A[95 - l] + it * 8);
					a = _mm256_mul_ps(a, a);
					__m256 r = _mm256_set1_ps(diagonal[l]);
					acc = _mm256_fmadd_ps(r, a, acc);
				}
				_mm256_storeu_ps(&m_invSymR[off + it * 8], acc);
			}
			off += 12 * 8;
			for (int it = 0; it < 12; it++)
			{
				int xBg = it * 8;
				for (int scan = xBg + 1; scan <= (96 - 8); scan++)
				{
					int lc = 96 - scan;
					__m256 acc = _mm256_setzero_ps();
					for (int l = 0; l < lc; l++)
					{
						__m256 a = _mm256_loadu_ps(A[95 - l] + xBg);
						__m256 b = _mm256_loadu_ps(A[95 - l] + scan);
						a = _mm256_mul_ps(a, b);
						__m256 r = _mm256_set1_ps(diagonal[l]);
						acc = _mm256_fmadd_ps(r, a, acc);
					}
					_mm256_storeu_ps(&m_invSymR[off], acc);
					off += 8;
				}
			}
#else
            //FIXME!
            float diagonal[96] = {};
            int off = block * triSz;
#endif
			for (int it = 0; it < 12; it++)
			{
				for (int lane = 0; lane < 7; lane++)
				{
					int xBg = it * 8 + lane;
					for (int h = (96 - 7 + lane); h < 96; h++)
					{
						int lc = 96 - h;
						float acc = 0.f;
						for (int l = 0; l < lc; l++)
						{
							float a = A[95 - l][xBg];
							float b = A[95 - l][h];
							float r = diagonal[l];
							acc += a * b * r;
						}
						m_invSymR[off] = acc;
						off++;
					}
				}
			}

		}
}

void SeSchwarzPreconditioner::BuildResidualHierarchy(const SeVec3fSimd* m_cgResidual)
{
	//printf("进入build阶段...\n");
	int nblock = (m_numVerts + 31) / 32; //level0簇的个数

	OMP_PARALLEL
	{
		std::vector<std::unordered_map<int,SeVec3fSimd>> diagTable(m_numLevel);
		
#pragma omp for

		for (int bid = 0; bid < nblock; ++bid)
		{
			for (int lane = 0; lane < 32; lane++)
			{
				int vid = lane + bid * 32;
				if (vid >= m_numVerts)
					break;
				//unsigned int vidOrig = m_MapperSortedGetOriginal[vid]; //原始索引位置
				int vidOrig = m_MapperSortedGetOriginal[vid];
				//printf("vidori %d\n", vidOrig);
				SeVec3fSimd r = m_cgResidual[vidOrig];  //结合原始索引从res列表中取值:长度应该=numvert 
				//printf("r is :%d\n", r.x);
				//std::cout << vid  << std::endl;
				//std::cout << vidOrig << std::endl;
				//std::cout << r.x << std::endl;
				//std::cout << r.y << std::endl;
				m_mappedR[vid] = r; //莫顿索引列表赋值
				if (m_numLevel > 1)
				{
					int next = m_goingNext[vid]; //level1的顶点索引
					//printf("next is : %d\n",next);
					//printf("now is : %f\n", m_mappedR[next].x);
					SeVec3fSimd temp = m_mappedR[next];
					//m_mappedR[next] = r + temp;
					m_mappedR[next] += r;  
					//printf("next now is : %f\n", m_mappedR[next].x);
				}
			}
		}
	}
		//我的理解：这里的mappedR就是右侧常数项映射到每个level所有顶点的结果,新level中顶点的residual项=上级level所有映射到该顶点的residual项的和
		if (m_numLevel > 2)
		{
			int bg = m_levelSize[1].y;
			int ed = m_levelSize[2].y;
			for (int vid = bg; vid < ed; vid++)
			{
				SeVec3fSimd r = m_mappedR[vid];
				int nx = vid;
				for (int lv = 2; lv < m_numLevel; lv++)
				{
					nx = m_goingNext[nx];
					m_mappedR[nx] += r;
				}
			}
		}
}

void SeSchwarzPreconditioner::SchwarzLocalXSym()
{
	constexpr int blockDim = 32;

	const int activeblockNum = m_totalNumberClusters / blockDim;

	const int triSz = (1 + 96) * 96 / 2 + 16 * 3;

	OMP_PARALLEL_FOR

		for (int blockIdx = 0; blockIdx < activeblockNum; ++blockIdx)
		{
			float cacheRhs[96];
			float cacheOut[96];
			for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + blockIdx * blockDim; //此处vid不再是vertex_index,是集成所有level的顶点索引{{0,1,2,3,4,5,6}{7,8}{9}}.{0,1,2,3,4,5,6}是原始顶点的索引{7,8}是level1的；
				auto rhs = m_mappedR[vid];   //叠加后的n个r:
				cacheRhs[laneIdx * 3 + 0] = rhs.x;
				cacheRhs[laneIdx * 3 + 1] = rhs.y;
				cacheRhs[laneIdx * 3 + 2] = rhs.z;
			}
#ifdef WIN32
			// ------------------------------------
			// diagnal
			int off = blockIdx * triSz;
			for (int it = 0; it < 12; it++)
			{
				__m256 diag = _mm256_loadu_ps(&m_invSymR[off + it * 8]); //invR是求得的逆矩阵
				__m256 diagRhs = _mm256_loadu_ps(&cacheRhs[it * 8]); //右侧常数项
				diag = _mm256_mul_ps(diag, diagRhs); //做矩阵乘法:M(-1)*b 应该得到精确解
				_mm256_storeu_ps(&cacheOut[it * 8], diag);
			}
			off += 8 * 12;
			for (int it = 0; it < 11; it++)
			{
				int xBg = it * 8;
				__m256 sigularRhs = _mm256_loadu_ps(&cacheRhs[xBg]);
				__m256 sigularResult = _mm256_setzero_ps();
				for (int scan = xBg + 1; scan <= (96 - 8); scan++)
				{
					__m256 mtx = _mm256_loadu_ps(&m_invSymR[off]);
					off += 8;

					__m256 scanRhs = _mm256_loadu_ps(&cacheRhs[scan]);
					sigularResult = _mm256_fmadd_ps(mtx, scanRhs, sigularResult);

					__m256 scanResult = _mm256_loadu_ps(&cacheOut[scan]);
					scanResult = _mm256_fmadd_ps(mtx, sigularRhs, scanResult);
					_mm256_storeu_ps(&cacheOut[scan], scanResult);
				}

				__m256 sigularResultOrg = _mm256_loadu_ps(&cacheOut[xBg]);
				sigularResult = _mm256_add_ps(sigularResultOrg, sigularResult);
				_mm256_storeu_ps(&cacheOut[xBg], sigularResult);
			}
#else
            //FIXME!
            int off = blockIdx * triSz;
#endif

			for (int it = 0; it < 12; it++)
			{
				for (int lane = 0; lane < 7; lane++)
				{
					int xBg = it * 8 + lane;
					float sigularRhs = cacheRhs[xBg];
					float singularResult = 0.0f;
					for (int h = (96 - 7 + lane); h < 96; h++)
					{
						float value = m_invSymR[off];
						off++;
						if (value)
						{
							singularResult += value * cacheRhs[h];
							cacheOut[h] += value * sigularRhs;
						}
					}
					cacheOut[xBg] += singularResult;
				}
			}


			for (int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + blockIdx * blockDim;

				SeVec3fSimd out(0.0f);
				out.x = cacheOut[laneIdx * 3 + 0];
				out.y = cacheOut[laneIdx * 3 + 1];
				out.z = cacheOut[laneIdx * 3 + 2];

				m_mappedZ[vid] = out;
				//std::cout << "out" << std::endl;
				//std::cout << out.x << std::endl;
			}

		}
}

void SeSchwarzPreconditioner::CollectFinalZ(SeVec3fSimd* m_cgZ)
{
	//如果不考虑碰撞之类的信息，单纯提取预处理器部分是否能有效？
	// 现在的工作重心可以简化到:给定一个矩阵，结合预处理器求其逆
	// 核心部分应该集中在buildlevelx,以及最后colectz，中间的碰撞可以暂时忽略
	OMP_PARALLEL_FOR
		//m_cgZ初始顶点位置信息
		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			//计算步骤: 1)莫顿索引vid->原始索引mapindex;2)根据vid索引得到映射信息z和table;3)将四层的map信息叠加;4)结果存储到m_cgz中
			unsigned int mappedIndex = m_MapperSortedGetOriginal[vid];//mappedIndex是最初始的顶点索引
			auto z = m_mappedZ[vid]; //初始值 z0 = r = b-A*0;z实质是更新的步长
			//SE::SeVec3fSimd z(0.0, 0.0, 0.0);
			//printf("%6f %6f %6f\n", z.x, z.y, z.z);
			//Int4 table = m_coarseTables[vid];//当前顶点的映射信息,index=6:tabel = [65,78,97],表示映射到不同层的索引
			//for (int l = 1; l < Math::Min(m_numLevel, 4); l++) 
			//{
			//	int now = table[l - 1];
			//	z += m_mappedZ[now]; 
			//}
			m_cgZ[mappedIndex] = z; //z表示步长不是更新后的顶点位置
		}
	//std::cout << "Compute finished" << std::endl;
}
