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
//�������ӦSeSchwarzPreconditionerͷ�ļ��е���&��������

#include "SeSchwarzPreconditioner.h" //ʩ�ߴ�Ԥ������
#include "SeIntrinsic.h" //Se intrinsic���ڵ�


SE_USING_NAMESPACE

void SeSchwarzPreconditioner::AllocatePrecoditioner(int numVerts, int numEdges, int numFaces)
{
	//Ԥ������䣺�����ڴ�֮���
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
	//const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets,  //��ײԪ���ж���������ṹ������
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

	//PrepareCollisionHessian(); //��ײ����
	PrepareHessian(diagonal, csrOffDiagonals, csrRanges); //�Խ�+ѹ���Խ�+ѹ����Χ
	LDLtInverse512(); //LDLt���棺����ǶԾ������������
}

void SeSchwarzPreconditioner::Preconditioning(SeVec3fSimd* z, const SeVec3fSimd* residual, int dim)
{
	//Ԥ�������ּ�
	Utility::MemsetZero(m_mappedZ); 
	Utility::MemsetZero(m_mappedR);

	BuildResidualHierarchy(residual); 

	SchwarzLocalXSym(); 
	 
	CollectFinalZ(z); 
}

void SeSchwarzPreconditioner::ComputeLevelNums(int bankSize)
{
	constexpr float SizeRatio = 1.5f; //banksize=32���ع�ģ

	//====��������������ʼ����

	int totalSz = 0;
	int nLevel = 1; //level1����

	//e.g������=64������õ��Ľ������2*32�����=65���������3*32������bank_size�Ĺ�ģ������ʽ
	int levelSz = (m_numVerts + bankSize - 1) / bankSize * bankSize; //�������������õ������ٳ���bank_size�õ�32�������Ĺ�ģ����Ϊԭʼ����������ܲ���32�ı����������޷����մؽ��л���
	totalSz += levelSz; //�ܹ��Ĺ�ģ=levelSz

	//�ײ����󣬵õ�levelSz = k*32,k��ʾ��ǰ��Ĵ��Ÿ���,32��Ӧ�صĹ�ģ;
	while (levelSz > 32) // 32��ʾ�����ڵ��ģ(һ�������ڵ�����Ľڵ���,��banksize=32��ͬһ����˼)
	{
		//node num of next level,nodesize=32,�ô�������һ����һ���ڵ��Ӧ�ϼ���һ����
		levelSz /= 32;
		nLevel++;
		levelSz = (levelSz + bankSize - 1) / bankSize * bankSize; //��ʱ�µ�level_size = 1*32<=32,ѭ����Ӧ����,3�������ڵ�ֻ�ܹ���һ���أ�����size=1*32;
		//Ҫע�����(levelSz + bankSize - 1) / bankSize����ʾ��ǰ��ڵ��ܹ��ɵĴ�������
		// (levelSz + bankSize - 1) / bankSize * bankSize - ��ʾ��ǰ���size����Ϊ������ܲ���32�ı����������size�Ƕ�ԭʼsize���е����Ľ�������漰��һ��level��������Ϣ��ֻ�ǶԵ�ǰ��ĵ���size
		totalSz += levelSz; 
	}

	m_numLevel = nLevel; //total level 
	m_totalSz = totalSz * SizeRatio; //Totalsize = totalsize*1.5?
}

void SeSchwarzPreconditioner::DoAlllocation()
{
	// do all location
	constexpr int bankSize = 32;	//==== size of diagonal block �Խǿ�Ĺ�ģ������32ָ���Ƕ���������Ǿ����еĸ���

	ComputeLevelNums(bankSize);

	//====һϵ�г�ʼ��resize����

	m_MapperSortedGetOriginal.resize(m_numVerts); //��������:��ԭʼ�ڵ�����������
	m_mapperOriginalGetSorted.resize(m_numVerts); //��������:Ī�ٿռ�����е�ԭʼ�ڵ����е���ӳ��
	m_mortonCode.resize(m_numVerts); //Ī�������͵��������洢Ī�ٱ��룬һ����ά�����ӦΨһ��Ī�ٱ���

	m_hessian32.Resize(bankSize, m_totalSz); //ָ�����У�32*(total_k*1.5*32), total_k��ԭʼ��������д�(����level0-3���е�)�ĸ���
	m_mappedZ.resize(m_totalSz); //map z total_size = total_k*32
	m_mappedR.resize(m_totalSz); //map r

	m_mappedPos.resize(m_numVerts); //����λ��ӳ��,��ģͬ����

	m_denseLevel.resize(m_numVerts); //densellevel
	m_goingNext.resize(m_numLevel * m_numVerts); //����*������
	m_coarseTables.resize(m_numVerts); //�ֲڱ���¼

	m_prefixOrignal.resize(m_numVerts);//ԭʼǰ׺
	m_fineConnectMask.resize(m_numVerts); //connect mask 
	m_nextConnectMsk.resize(m_numVerts); 
	m_nextPrefix.resize(m_numVerts);


	int triSz = (1 + 96) * 96 / 2 + 16 * 3; //�洢���������4704
	const int blockNum = m_totalSz / 32; //block_num�� ��1.5total_size�������blocknum = ԭʼnum*1.5
	m_invSymR.resize(triSz * blockNum); //ϵͳ������棺
	m_additionalHessian32.resize(m_totalSz + 1);
	//����hessian��,��С=totalSz+1,32ֻ�Ƕ����ģ,��ÿ������������3�����깹�ɵ�
	//m_additionalHessian32����totalSz+1��3*3�Խ��󹹳ɵ�����

	//====

	m_levelSize.resize(m_numLevel + 1); //��level+1,�ӵ�0��,�洢��ԭʼ������
	m_CoarseSpaceTables.Resize(m_numLevel, m_numVerts); //�ֲڿռ����ͬlevel�Ͷ��������0����Ӧ��ʼ״̬

	int maxNeighbours = 0;

	for (int vId = 0; vId < m_numVerts; ++vId)
	{
		maxNeighbours = Math::Max(m_neighbours->Size(vId) + 1, maxNeighbours);//�������нڵ㣬�õ���������������ڵĶ������+1(��������)
	}

	//һϵ��resize����
	m_mappedNeighborsNum.resize(m_numVerts);
	m_mappedNeighborsNumRemain.resize(m_numVerts);//remain

	//������Сдresize,����������Resize
	m_mappedNeighbors.Resize(maxNeighbours, m_numVerts); //�洢���������ھ�����:[i,j]�洢��j������ĵ�i���ھ��±�
	m_mappedNeighborsRemain.Resize(maxNeighbours, m_numVerts);

	const int maxCollisionPerVert = 32; //�����ײ������ֻ����������ײ
	m_maxStencilNum = m_numVerts * maxCollisionPerVert;
	m_stencils.resize(m_maxStencilNum);
	m_stencilNum = 0;
}

void SeSchwarzPreconditioner::ComputeTotalAABB()
{
	m_aabb = SeAabb<SeVec3fSimd>(); //��ʼ��һ����Χ��

	ComputeAABB(); // TODO: omp reduce
	//printf("BVH total end!\n");
}


void SeSchwarzPreconditioner::ComputeAABB()//������Χ����� axis alien bound box
{

//MSVC compiler is stuck with OpenMP version 2.0, and unfortunately for you, reduction(max:) was only introduced with version 3.1 of the OpenMP C/C++ standard (that was in September 2011) 

//#pragma omp parallel for reduction(aabbAdd:m_aabb)
	//printf("BVH action!\n");
	for (int vid = 0; vid < m_numVerts; ++vid) // need omp custom reduction operation
	{
		//printf("inital lower : %f\n", m_aabb.Upper.x);
		m_aabb += m_positions[vid];//�������㣬�õ��������ж����һ����Χ��

	}
	//printf("BVH end!\n");
}

void SeSchwarzPreconditioner::SpaceSort()
{
	FillSortingData(); //�ռ����򣺽��Ī����Ͷ�����������
	DoingSort();
}

void SeSchwarzPreconditioner::FillSortingData()
{
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			SeVec3fSimd temp = (m_positions[vid] - m_aabb.Lower) / m_aabb.Extent(); //��ʱ����������������������[0,1]֮��

			SeMorton64 mt;

			mt.Encode(temp.x, temp.y, temp.z); //�õ���ǰ�����Ī����

			m_mortonCode[vid] = mt; //Ī�������洢

			m_MapperSortedGetOriginal[vid] = vid; //ԭʼ��������Ƕ�������������֮����Ҫ��������ʼ[0,1,2,3,4,5,6]
		}
}


void SeSchwarzPreconditioner::DoingSort()
{
	//ͨ��Ī����õ�ԭʼ��������index��Ī�ٶ����value��ԭʼ����
	//==== sort m_MapperSortedGetOriginal by m_mortonCode1; ͨ��Ī���뽫ԭʼ������������(��С����)
	//[0,1,2,3,4,5,6] -> [1,3,0,6,5,4,2]����ʾ����i�����Ľ����[1,3,0,6,5,4,2]; ��value,index��ԭʼ��value��Ԫ����Ī����ĵ�index��Ԫ��
	//ʵ�ʲ����������������б�[1,3,0,6...]�����ҵ�ԭʼ���Ϊ1,3,0,6�Ķ��㰴˳�����±���Ϊ0,1,2...
	std::sort(m_MapperSortedGetOriginal.begin(), m_MapperSortedGetOriginal.end(), [&](int i, int j) { return m_mortonCode[i] < m_mortonCode[j]; }); //ǰ��������ָ����Χ��������ָ������Ĺ���
}

void SeSchwarzPreconditioner::ComputeInverseMapper()
{
	//��index��ԭʼ�����value��Ī������������
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			int originalIndex = m_MapperSortedGetOriginal[vid];

			m_mapperOriginalGetSorted[originalIndex] = vid; //���ϱ߶�Ӧ�ľ���[2,0,6,1,5,4,3],��ʾĪ�����еĵ�value��Ԫ����ԭʼ��index��Ԫ�أ�value,index����ӳ��
		}
}


void SeSchwarzPreconditioner::MapHessianTable()
{
	OMP_PARALLEL_FOR

		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			const auto  originalIndex = m_MapperSortedGetOriginal[vid]; //�õ�ԭʼ��Ī�ٹ����С�����index

			m_mappedPos[vid] = m_positions[originalIndex]; //���������긳ֵ��Ī�ٶ���λ��

			//====

			int numNeighbors = m_neighbours->Size(originalIndex); //�ھӸ�����ֵ

			const int* neighbors = m_neighbours->IdxPtr(originalIndex); //��ǰ���������ھӵ��������б�

			m_mappedNeighborsNum[vid] = numNeighbors + 1;	// include self ���б�洢����������ھӶ������

			m_mappedNeighbors[0][vid] = vid;				// include self �洢����������ھӶ��������,��vid������ĵ�0���ھ����Լ�

			for (int k = 1; k < numNeighbors + 1; ++k)
			{
				unsigned int originalNeighbors = neighbors[k - 1];//m_neighbors[k][originalIndex];

				m_mappedNeighbors[k][vid] = m_mapperOriginalGetSorted[originalNeighbors];
			}
		}
}

void SeSchwarzPreconditioner::MapCollisionStencilIndices()
{
	//��ײ�����Ī�ٿռ��ӳ��
	m_stencilIndexMapped.resize(m_stencilNum);

	OMP_PARALLEL_FOR

		for (int i = 0; i < m_stencilNum; i++)
		{
			const auto& s = m_stencils[i];

			for (int vi = 0; vi < s.verextNumPerStencil; vi++)
			{
				//��Ī�������¶���ײ����:��ײ��1(vf�ṹ)��ԭʼindex��[1,3,5,6] -> ��Ī�ٿռ���� [4,7,2,8]
				//m_stencilIndexMapped[i]��ʾ��Ī�ٿռ��е�i����ײ�����index�ṹ
				m_stencilIndexMapped[i][vi] = m_mapperOriginalGetSorted[s.index[vi]]; //ת��ӳ�� 
			}
		}
}

void SeSchwarzPreconditioner::PrepareCollisionStencils(const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts)
{
	//�ο����ף�PSCC:������֮���ee&vfһ����15������
	/*ʵ��:4������, 5����, 2����ļ򵥰���,{0,1,2,3};efset=[{00}{10}{20}{21}{31}{41}];eeset=[{01}{02}{03}{12}{14}{23}{24}{34}];vfset=[{00}{10}{20}{11}{21}{31}],efcount=[],eecount=[],vfcount=[]*/

	//������ײ������\����\���� 
	int efNum = efCounts[m_numEdges]; //count����,�ۼ��б�,num��Ӧ����
	int eeNum = eeCounts[m_numEdges];
	int vfNum = vfCounts[m_numVerts];

	//���п�����ײ������
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

			//������ײ������б�: ����,һ������;��,��������(�˵�);��,��������(����);
			//e.g:{ee:ÿ������Ҫ������������ȷ��,����verextNumPerStencil=2+2}��{ef:,����verextNumPerStencil=2+3}
			if (i < efNum) //process ef
			{
				const auto& pair = efSets[i]; //��i��e-f���Ӵ�

				if (pair.m_eId < 0 || pair.m_fId < 0) { continue; }

				const auto& edge = m_edges[pair.m_eId];
				const auto& face = m_faces[pair.m_fId];

				s.verextNumPerStencil = 5;
				s.vertexNumOfFirstPrimitive = 2;

				s.index[0] = edge[0]; //edge:�������ݣ�һ�����������˵㹹��
				s.index[1] = edge[1];
				s.index[2] = face[0]; //������
				s.index[3] = face[1];
				s.index[4] = face[2];

				//Ȩ��ֵ�����ĵ������й�
				s.weight[0] = pair.m_bary[0]; //bary����,�������� bary.x,.y,.z
				s.weight[1] = 1.f - pair.m_bary[0];
				s.weight[2] = -pair.m_bary[1];
				s.weight[3] = -pair.m_bary[2];
				s.weight[4] = -(1.f - pair.m_bary[1] - pair.m_bary[2]);

				s.direction = pair.m_normal; //�����뷨��һ��

				s.stiff = pair.stiff; //�ն�ϵ��

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
				s.vertexNumOfFirstPrimitive = 3; //ԭʼ����num:��-����,��Ϊԭʼ = 3
				
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

			int stencilNum = Intrinsic::AtomicAdd(&m_stencilNum, 1); // �ȼ���ִ��stencilNum = m_stencilNum;m_stencilNum+=1;

			m_stencils[stencilNum] = s; //����ײ�ṹ�洢��m_stencils��
		}

	MapCollisionStencilIndices(); //1
}

void SeSchwarzPreconditioner::ReorderRealtime()
{
	//ʵʱ����
	Utility::MemsetZero(m_levelSize);

	BuildConnectMaskL0(); //��level-0��ʼ

	BuildCollisionConnection(m_fineConnectMask.data(), nullptr); //2

	PreparePrefixSumL0(); 

	BuildLevel1(); 

	for (int level = 1; level < m_numLevel; level++)
	{
		//����level
		Utility::MemsetZero(m_nextConnectMsk);

		BuildConnectMaskLx(level); //����mask

		BuildCollisionConnection(m_nextConnectMsk.data(), m_CoarseSpaceTables[level - 1]); //��ײ��ϵ

		NextLevelCluster(level); //

		PrefixSumLx(level); //

		ComputeNextLevel(level); //������һ��level
	}

	TotalNodes(); //�ڵ�����

	AggregationKernel(); //���ٺ�
}

void SeSchwarzPreconditioner::BuildConnectMaskL0()
{
	//������ȷ������connetmask�Ǻκ���? ָ����������֮������ӹ�ϵ:��ͬ���ڵ�������Ϣ�Ͳ�ͬ���ڵ�������Ϣ(����level����,��Ҫ���Ƕ���ϲ���Ϣ���µ�ǰlevel�Ķ���������Ϣ)
	//mask��¼�˴��ڶ���֮���������Ϣ�������Ǵ��⣡

	//��ÿ��ԭʼ������level0�ִأ���������mask��neighbor��Ϣ
	//L0����ϵmask? �ⲿ��Ӧ�ö�Ӧ��ԭʼ��
	const unsigned int warpDim = 32; //Ť��ά�� = banksize
	int nWarps = (m_numVerts + warpDim - 1) / warpDim; //������ warp:��,��������Ϊ������

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < nWarps; ++warpIdx) 
		{	// warp-id:�ڼ����� 0-nwarp
			for (int laneIdx = 0; laneIdx < warpDim; ++laneIdx)
			{
				//laneid:���ڵڼ�������0-31
				int vId = warpIdx * warpDim + laneIdx; // vID: 0-nWarps*warpdim/numvert

				if (vId >= m_numVerts)
				{
					//����������������ѭ�����������������32�ı����п��ܷ���
					break;
				}
				int numNeighbors = m_mappedNeighborsNumRemain[vId]; //��ǰ������ھӸ���,mapped˵����Ī�ٿռ��ϲ���
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
				2147483648 2^0 - 2^31 �ڲ�ѭ��һ�εõ���ʼ��mask��ʽ
				*/
				unsigned int connetMsk = 1U << laneIdx;	// include self 1U-32λ�޷�������,����idx��������λ

				int nk = 0;

				for (int k = 0; k < numNeighbors; k++)
				{
					int vIdConnected = m_mappedNeighborsRemain[k][vId]; //�õ���k���ھӵ�����

					int warpIdxConnected = vIdConnected / warpDim;//�õ��ص���������ʾvID����������ӵ���Ϣ:e.g vID=1,������34�Ķ�����������vID���1������

					if (warpIdx == warpIdxConnected) // in the same warp ����ǰ���������ھ���ͬһ������
					{
						unsigned int laneIdxConnected = vIdConnected % warpDim; //ȡ����

						//����mask,��λ���ֵ,������Ӧ�õ���0
						connetMsk |= (1U << laneIdxConnected);  //��ͬ���������ھ�:mask = (1U << laneIdx) | (1U << laneIdxConnected),���:��һ����1��λ=1;����0��λ=0,��ֵ����
						//�پٸ�����:�Ա��Ϊ0���ʵ�Ϊ������ʼmask=1U<<0 = 0001 = 1;�ش�С��4(������λֻ��Ҫ4λ����)�����ڵ�����{0,2,6},{0,2}��0����ͬ��,����ѭ��,0%4=0,2%4=2,����mask = 0001|0100 = 0101 
						//�۲춥���connetmask����ֱ�ӹ۲쵽����(ֻ���Ǵ���)�����������Ϣ
					}
					else
					{
						m_mappedNeighborsRemain[nk][vId] = vIdConnected; //��ǰ����������ӵ��Խ��ͬ�أ�����Ϊremain�ھ�(���µ��ھ�),�����������´��浽remainneighbor��
						nk++;
					}
				}
				m_mappedNeighborsNumRemain[vId] = nk; //���µ�ǰ�ڵ���ھ���Ϣ����֮ǰ��������,���ٸ��� = �뵱ǰ������ͬһ�����������ӵĶ���

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
				
				//finemask��¼��ÿ�����������ڴ����������������Ϣ��������
				//m_mappedNeighborsRemain[:][vertexindex]ÿ���������������ж����������Ϣ:remainneighbor����������
				//�����ںʹ���������������ֲ�ͬ������������� - ���������������=ԭʼ�ھӸ���
				m_fineConnectMask[vId] = connetMsk; //�洢ÿ���ڵ��mask��Ϣ,fine�洢���ȫ����mask���,����ͬVertnum
			}
		}
}

void SeSchwarzPreconditioner::BuildCollisionConnection(unsigned int* pConnect, const int* pCoarseTable)
{
	//������ײ��ϵ:��ײ������ʱ����
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
					if (myID == otID) // redundant ? ����ģ�
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
	//�˴�ʵ�������Ӷ���֮�����ת:�ֱ�ǰ�ڵ���ת���������ӵ�ÿһ�����㣬������mask,��֤���ӵĶ��㹹��һ���ջ�
	
	//����: �ع�ģ(dim),������
	int warpDim = 32;
	int nWarps = (m_numVerts + warpDim - 1) / warpDim;

	OMP_PARALLEL_FOR

		for (int warpIdx = 0; warpIdx < nWarps; ++warpIdx)
		{
			unsigned int cacheMask[32]; //��ʱ�������:�洢mask��Ϣ;����32,�ⲿѭ��һ�θ���һ��

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vId = laneIdx + warpIdx * warpDim; //��ǰ��������

				if (vId >= m_numVerts) { break; }

				cacheMask[laneIdx] = m_fineConnectMask[vId]; //�������ݵ�cache or cache�б�ֵ
			}

			Intrinsic::SyncWarp(); //CUDA�߳�ͬ��

			int prefixSum = 0; //��������

			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * warpDim; //��ʵ���������

				if (vid >= m_numVerts) { break; }

				unsigned int connetMsk = cacheMask[laneIdx];//��ȡcache�����ݡ���connet

				unsigned int visited = (1U << laneIdx); //��λ����: ��1U����landIdxλ,����Ľ��=1U*(2^(laneIdx))

				//unsigned int = -1,��ʾ�޷������͵����ֵ = 11111111111111111111111111111111,��ʱ˵����ǰ������������ж�������
				
				//��ѭ��Ϊ�˱�֤�����ڵ������ӣ��������ӵĶ���֮�䱣֤��ͨ������mask�ٴθ���
				while (connetMsk != -1) 
				{
					//����connetmsk�з�0��ÿһλ,visited��¼����״̬,todoָʾ��������ʾ�Ƿ������ϣ�
					unsigned int todo = visited ^ connetMsk; // �Զ����Ƶ�ÿһλִ����ͬλ=0;��ͬλ=1

					if (!todo)
					{
						//todo = 0����ѭ��:��visited == connetMsk ÿһλ�����
						break;
					}

					unsigned int nextVisit = Intrinsic::Ffs(todo) - 1; //����todo�׸�����λ������

					visited |= (1U << nextVisit); //��ת����һ�����ӵĵ㣬����visit״̬

					connetMsk |= cacheMask[nextVisit]; //����mask,Ϊ�˱�֤�����������ӵ�����ͨ��
				}
				m_fineConnectMask[vid] = connetMsk; //�洢���º��mask��Ϣ

				//__popcount(temp)����temp������λ��=1�ĸ���; &��0��0;|��1��1(����Ҳ��1)
				//LanemaskLt:(1U << laneId) - 1 = 2^(laneId)-1;����:����connetMskǰlaneIdx������λ��=1�ĸ���
				//�˴�ִ�еĹ��ܴ��ɣ����㷵�صĽ��ǡ���Ƕ�����λ����λ����+1,�п���ֱ�ӵ��õĺ���Ϊʲô��Ҫ�������
				unsigned int electedPrefix = Intrinsic::Popcount32(connetMsk & Intrinsic::LanemaskLt(laneIdx)); 
				if (electedPrefix == 0)
				{
					//=0��ʾ��ǰid�����㶼������,sum+1
					prefixSum++;
				}
			}
			m_prefixOrignal[warpIdx] = prefixSum; //�׸����ӵĶ����λ�ã�
		}
}

void SeSchwarzPreconditioner::BuildLevel1()
{
	constexpr unsigned int blockDim = 1024;//32*32

	int nBlocks = (m_numVerts + blockDim - 1) / blockDim;//level1������ = level0/32��level1�Ĵ�����= level0/32/32

	OMP_PARALLEL_FOR
	// �ṹ:for{ for1{ for11};for2{};for3{for31 for32 for33 for34}   }

		//����level1�Ĵ�:block>(super node)>(node)
		for (int blockIdx = 0; blockIdx < nBlocks; ++blockIdx) //������
		{
			unsigned int warpPrefix[32];//��ʱ���������ڱ���ÿ�����ڵ�prefix��Ϣ,����32,ÿ�δ�ѭ������

			unsigned int theBlockGlobalPrefix = 0; //blockȫ��ǰ׺

			unsigned int globalWarpOffset = blockIdx * 32; //��ƫ��:һ����ƫ��32������

			//������ȷ���庬��,�õ�theBlockGlobalPrefix
			for (int warpIdx = 0; warpIdx < blockIdx; ++warpIdx) //
			{
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;
					theBlockGlobalPrefix += m_prefixOrignal[threadIdx]; //����,���������µ��ӵ���Ϊ����������ǳ���
				}
			}

			//��prefix��Ϣ��������ʱ�б���
			for (int warpIdx = 0; warpIdx < 32; ++warpIdx)
			{
				warpPrefix[warpIdx] = m_prefixOrignal[globalWarpOffset + warpIdx];
			}

			Intrinsic::SyncBlock();

			//����coarseTable & goingMap
			for (unsigned int warpIdx = 0; warpIdx < 32; ++warpIdx) //32-���ڱ���
			{
				unsigned int theWarpGlobalPrefix = theBlockGlobalPrefix;

				for (unsigned int prevWarpIdx = 0; prevWarpIdx < warpIdx; ++prevWarpIdx)
				{
					theWarpGlobalPrefix += warpPrefix[prevWarpIdx]; //�ٴε���
				}

				//����electedmask
				unsigned int electedMask = 0;
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32; //level1�³�������ı��(thread-supernode)

					int vid = threadIdx + blockIdx * blockDim; //initial vertex index

					if (vid >= m_numVerts) { break; }


					if (vid == m_numVerts - 1) //���һ������
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

				//����lanemake[],����32
				unsigned int lanePrefix[32];
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= m_numVerts) { break; }


					lanePrefix[laneIdx] = Intrinsic::Popcount32(electedMask & Intrinsic::LanemaskLt(laneIdx));

					lanePrefix[laneIdx] += theWarpGlobalPrefix;
				}

				//����coarseTable & goingnext
				for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
				{
					unsigned int threadIdx = laneIdx + warpIdx * 32;

					int vid = threadIdx + blockIdx * blockDim;

					if (vid >= m_numVerts) { break; }


					unsigned int connMsk = m_fineConnectMask[vid];

					unsigned int elected_lane = Intrinsic::Ffs(connMsk) - 1;

					unsigned int theLanePrefix = lanePrefix[elected_lane];

					m_CoarseSpaceTables[0][vid] = theLanePrefix;

					m_goingNext[vid] = theLanePrefix + (m_numVerts + 31) / 32 * 32;  //goingNext������ÿһ�����ʵ㵽��һ����ӳ��(ע���ǰ�������level�����ж��㣡)
				}
			}
		}
}


void SeSchwarzPreconditioner::BuildConnectMaskLx(int level)
{
	//������ͬ��֮�����ϵ
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
	//��һ���Ĵؼ���
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
	// ǰ׺���
	const int levelNum = m_levelSize[level].x;

	const int levelBegin = m_levelSize[level].y;

	unsigned int blockDim = 1024;//32*32

	int nBlocks = (m_numVerts + blockDim - 1) / blockDim; //�������һ��Ĵصĸ���

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
		//�����¸�level
		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			int next = m_CoarseSpaceTables[level - 1][vid];

			m_CoarseSpaceTables[level][vid] = m_nextConnectMsk[next];
		}
}

void SeSchwarzPreconditioner::TotalNodes()
{
	//��i����ܵĴظ���
	m_totalNumberClusters = m_levelSize[m_numLevel].y;
	//printf("%d\n", m_totalNumberClusters);
}

void SeSchwarzPreconditioner::AggregationKernel()
{
	const int warpDim = 32;

	int numWarps = (m_numVerts + warpDim - 1) / warpDim; //��һ����cluster����

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

				curr[laneIdx] = vid; //����32,�����ڴ��ڵı���

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
					//e.g:�����Ƕ���0,curr[0]=32��ʾ��ǰ���ڵ�0�������ǳ�ʼ��32������;m_goingNext[curr[0]]=66��next[0]=1��ʾ��ǰ���ڵ�0����������һ��level����Ϊ7�Ķ���
					//�൱��:��laneidxΪ����,curr�洢���ڵ�ǰlevel��Ӧ�Ķ�������;next�洢��ӳ�䵽��һ���Ķ�������
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
					coarseTable[laneIdx][l] = next[laneIdx]; //������L����1
				}
			}


			for (unsigned int laneIdx = 0; laneIdx < 32; ++laneIdx)
			{
				int vid = laneIdx + warpIdx * 32;

				if (vid >= m_numVerts) { break; }

				m_denseLevel[vid] = aggLevel[laneIdx]; //agglevel&denselevel��֮�������ò�ƶ�û���õ�;�ݲ�����

				m_coarseTables[vid] = coarseTable[laneIdx]; 
				//coarse����Ϊ32;m_coarseTables����Ϊnumvert:��¼��ÿ���������е�ӳ����Ϣ
				//���룺��ʽ����m_coarseTable = [[65,66,67],[65,66,67],[65,66,67]]...��������Ϊ1���б�[0,0,0],��ʾ��ʼ��һ������ֱ�ӳ�䵽level1,2,3��������65,66,67;
			}
		}
}

void SeSchwarzPreconditioner::AdditionalSchwarzHessian2(SeMatrix3f hessian, std::vector<SeMatrix3f>& pAdditionalHessian, SeArray2D<SeMatrix3f>& pDenseHessian, int v1, int v2, const std::vector<int>& pGoingNext, int nLevel, int vNum)
{
	//�����������ã�����ͬһ���ڽڵ����ײ
	//1) additional ����������ڶ����ܹ��ϲ�������levelӳ���ĳ������㣬��ӳ�䵽��һ��level�ĳ�������(�Ա���һ����ӳ����һ��)��Ȼ����и�ֵ�������ⲿ�ֲ��������һ�㸳ֵ����Ϊ���һ���޷���ӳ�䵽��һ�㣡
	//2) dense ԭʼ��ײ�������ڵ�һֱѭ���������ܹ��ϲ��������ڶ����ܹ��ϲ�Ϊһ��������level��ӳ���ĳ����ڵ�����������hessian32���и�ֵ(��¼��ԭʼ�ڵ�ӳ�䵽ĳ��levelͬһ�����ڳ����ڵ����ײ)������һ������



	int level = 0;
	unsigned int myID = v1;
	unsigned int otID = v2;
	const unsigned int bank = 32; //�صĴ�С

	while (myID / bank != otID / bank && level < nLevel)
	{
		//�����㲻�ܻ���ͬһ����������һ���ڵ�ֱ�����߿��Ժϲ���ͬһ������ߵ����������
		myID = pGoingNext[myID];
		otID = pGoingNext[otID];
		level++;
	}

	if (level >= nLevel) //����ڶ���ѭ����ֹ������return
		return;
	
	//����level < nLevel && myID / bank == otID / bank
	//��m_hessian32��3*3�Խǿ�������: ÿ��level�����ܸ���0-3
	Intrinsic::AtomicAdd(&pDenseHessian[otID % bank][myID], hessian); //hessian32���յ�hessian�����ڶ��ߵ�ǰ�ܹ��ϲ�����ֵ
	Intrinsic::AtomicAdd(&pDenseHessian[myID % bank][otID], hessian);

	if (level < nLevel - 1)
	{
		//��level1-3����
		myID = pGoingNext[myID];
		otID = pGoingNext[otID];

		if (myID == otID) //������������һ��level�ܹ��ϲ�Ϊһ���ڵ㣬����64��65�ϲ�Ϊ100
		{
			Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian * 2.0f); //����ײ�����ܹ��ϲ��ĳ������㸳ֵ2*hessian
		}
		else
		{
			//����vnumҲû����
			//�˴�������ô���ܲ�����أ��벻ͨ = =�������Ҳ���ԣ�banksize != supernode size
			//����˵�����else����̫�ԣ�Ӧ�ö�Ӧ�����Ǹ�if?
			std::cout << "Label this running!" << std::endl;
			Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian);
			Intrinsic::AtomicAdd(&pAdditionalHessian[otID], hessian);
		}
	}
	//������������ڣ�level=3<nlevel=4ʱ������ϱ�ѭ������ myID / bank == otID / bank,˵�����������һ��պ���һ�����ڣ������һ�㲻�ٺϲ����㣻�ǰ����������ͻ�ȷʵ���һ�㳬������ĸ�ֵ��
	//else 
	//{
	// //�������Ҳ���ģ���ע�ͣ�elseҲ����level = nLevel - 1,�������һ��û��goingnext���������ڵ���һ�����ڣ�������ͬһ���ڵ㣬��ʱ����Ҫ���ж��⸳ֵ
	//	Intrinsic::AtomicAdd(&pAdditionalHessian[myID], hessian);
	//	Intrinsic::AtomicAdd(&pAdditionalHessian[otID], hessian);
	//}
	//20220915:˼������ԭʼ��������ȷ�ģ����Ǵ˴�ΪʲôҪ���ӣ����ӵ����������������ڵ�ǰ��λ��ͬһ���أ��������������Ӧ����һ������������������ֱ𸽼�һ��hessian���󣻵����������Ѿ�λ�����һ��ʱ������������һ��������Ҳ��������һ���ĸ���
	//���ڵ�ǰ����ĸ����Ѿ���ֵ����pdensehessian�����У��������ⲻ��


}

void SeSchwarzPreconditioner::PrepareCollisionHessian()
{
	//ִ����ײ��Ӧ���ڳ�ʼϵͳ�������޸�,�˴�����ײ��Ӧ���ڼ�⵽��ײʱ������һ��ʱ������ʩ��һ��������ֹ�����ײ�ķ���,����Ϊ������ԭʼϵͳ������
	//����ų��������Ӵ���ײ����wang2020 repulsion
	OMP_PARALLEL_FOR
		//m_additionalHessian32��ʼ��=0
		for (int i = 0; i < m_stencilNum; i++)
		{
			const auto& s = m_stencils[i];
			auto idx = m_stencilIndexMapped[i]; //ʹ��map���idx

			Float3 d(s.direction[0], s.direction[1], s.direction[2]); //��s�ķ���ֵ��d

			SeMatrix3f hessian = OuterProduct(d, d * s.stiff); //�õ�3*3���� �� hessian,��������ϵ��


			for (int it = 0; it < s.verextNumPerStencil; it++)
			{
				
				Intrinsic::AtomicAdd(&m_additionalHessian32[idx[it]], hessian * Math::Square(s.weight[it])); //���θ���m_additionalHessian32[idx[it]]����Ī�ٿռ����������
			}
			
			for (int ita = 0; ita < s.verextNumPerStencil; ita++)
			{
				for (int itb = ita + 1; itb < s.verextNumPerStencil; itb++)
				{
					//���θ���
					AdditionalSchwarzHessian2(s.weight[ita] * s.weight[itb] * hessian, m_additionalHessian32, m_hessian32, idx[ita], idx[itb], m_goingNext, m_numLevel, m_numVerts);
				}
			}
		}
}

void SeSchwarzPreconditioner::PrepareHessian(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges)
{

	const unsigned int bank = 32;

	int nVC = (m_numVerts + 31) / 32 * 32; //��level1������

	OMP_PARALLEL_FOR

		for (int vid = nVC; vid < m_totalNumberClusters; ++vid) //m_totalNumberClusters���нڵ�ĸ���: ���п��ܰ��������ڵĽڵ㣺��Ϊ����32�ı��ȴճ�32
		{
			auto oldDiagonal = m_additionalHessian32[vid]; //3*3
			int myID = vid;
			Intrinsic::AtomicAdd(&m_hessian32[myID % bank][myID], oldDiagonal); //���ӵ�ԭʼ�Խ���
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
		//���̴߳���
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
				int oldNum = m_mappedNeighborsNum[vid]; //���������neighbor����

				auto oldDiagonal = diagonal[vidOrig] + m_additionalHessian32[vid];
				m_hessian32[vid % bank][vid] += oldDiagonal;   // self diagonal

				for (int k = 1; k < oldNum; k++)  //index��1��ʼ,�������Լ�����ھ�
				{
					const unsigned int neighbor = m_mappedNeighbors[k][vid];
					SeMatrix3f mat = csrOffDiagonals[csrRanges[vidOrig] + k - 1]; //csrrange����start,�洢ÿ���ʵ��׸�offdiagnoal���б�csrrange[]�е�index = ǰi-1��������ھ�֮��(�����Լ�)
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
					if (level <= 1) // ����block�ֲ��У�level 1��λ��Ҳһ�����ڱ��߳�
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
	//LDlt(LDLת��)����
	const int triSz = (1 + 96) * 96 / 2 + 16 * 3;

	const int activeblockNum = m_totalNumberClusters / 32; //m_totalNumberClustersȫ���ڵ�ĸ���(��������ڵ�)��activeblockNum������ĸ���

	OMP_PARALLEL_FOR
		//printf("%d\n",CPU_THREAD_NUM);
		for (int block = 0; block < activeblockNum; block++)
		{
			float A[96][96] = {};
			for (int x = 0; x < 32; x++)
			{
				for (int y = 0; y < 32; y++)
				{
					SeMatrix3f temp = m_hessian32[y][x + block * 32]; //m_hessian = (����ײϵ����+������ײ����)+�༶�򹹽��Ľ�� = final ϵ������ = �����õ���ĸ���(totalclusters)��32*32����; 

					if (x == y && temp(0, 0) == 0.0f)
					{
						temp = SeMatrix3f::Identity(); //��λ��
					}
					for (int ii = 0; ii < 3; ii++)
					{
						for (int jj = 0; jj < 3; jj++)
						{
							A[x * 3 + ii][y * 3 + jj] = temp(ii, jj); //��A������,�Ҳ���ӵ�λ����
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

			// ������Ԫ
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
	//printf("����build�׶�...\n");
	int nblock = (m_numVerts + 31) / 32; //level0�صĸ���

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
				//unsigned int vidOrig = m_MapperSortedGetOriginal[vid]; //ԭʼ����λ��
				int vidOrig = m_MapperSortedGetOriginal[vid];
				//printf("vidori %d\n", vidOrig);
				SeVec3fSimd r = m_cgResidual[vidOrig];  //���ԭʼ������res�б���ȡֵ:����Ӧ��=numvert 
				//printf("r is :%d\n", r.x);
				//std::cout << vid  << std::endl;
				//std::cout << vidOrig << std::endl;
				//std::cout << r.x << std::endl;
				//std::cout << r.y << std::endl;
				m_mappedR[vid] = r; //Ī�������б�ֵ
				if (m_numLevel > 1)
				{
					int next = m_goingNext[vid]; //level1�Ķ�������
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
		//�ҵ���⣺�����mappedR�����Ҳೣ����ӳ�䵽ÿ��level���ж���Ľ��,��level�ж����residual��=�ϼ�level����ӳ�䵽�ö����residual��ĺ�
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
				int vid = laneIdx + blockIdx * blockDim; //�˴�vid������vertex_index,�Ǽ�������level�Ķ�������{{0,1,2,3,4,5,6}{7,8}{9}}.{0,1,2,3,4,5,6}��ԭʼ���������{7,8}��level1�ģ�
				auto rhs = m_mappedR[vid];   //���Ӻ��n��r:
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
				__m256 diag = _mm256_loadu_ps(&m_invSymR[off + it * 8]); //invR����õ������
				__m256 diagRhs = _mm256_loadu_ps(&cacheRhs[it * 8]); //�Ҳೣ����
				diag = _mm256_mul_ps(diag, diagRhs); //������˷�:M(-1)*b Ӧ�õõ���ȷ��
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
	//�����������ײ֮�����Ϣ��������ȡԤ�����������Ƿ�����Ч��
	// ���ڵĹ������Ŀ��Լ򻯵�:����һ�����󣬽��Ԥ������������
	// ���Ĳ���Ӧ�ü�����buildlevelx,�Լ����colectz���м����ײ������ʱ����
	OMP_PARALLEL_FOR
		//m_cgZ��ʼ����λ����Ϣ
		for (int vid = 0; vid < m_numVerts; ++vid)
		{
			//���㲽��: 1)Ī������vid->ԭʼ����mapindex;2)����vid�����õ�ӳ����Ϣz��table;3)���Ĳ��map��Ϣ����;4)����洢��m_cgz��
			unsigned int mappedIndex = m_MapperSortedGetOriginal[vid];//mappedIndex�����ʼ�Ķ�������
			auto z = m_mappedZ[vid]; //��ʼֵ z0 = r = b-A*0;zʵ���Ǹ��µĲ���
			//SE::SeVec3fSimd z(0.0, 0.0, 0.0);
			//printf("%6f %6f %6f\n", z.x, z.y, z.z);
			//Int4 table = m_coarseTables[vid];//��ǰ�����ӳ����Ϣ,index=6:tabel = [65,78,97],��ʾӳ�䵽��ͬ�������
			//for (int l = 1; l < Math::Min(m_numLevel, 4); l++) 
			//{
			//	int now = table[l - 1];
			//	z += m_mappedZ[now]; 
			//}
			m_cgZ[mappedIndex] = z; //z��ʾ�������Ǹ��º�Ķ���λ��
		}
	//std::cout << "Compute finished" << std::endl;
}
