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
#include "SeMorton.h" //Ī����
#include "SeAabbSimd.h" //aabbximd
#include "SeCollisionElements.h" //��ײ


SE_NAMESPACE_BEGIN

class SeSchwarzPreconditioner
{

public:
	
	//==== input data
	//��ʼָ��һ����ָ��
	const SeVec3fSimd* m_positions;		// the mesh nodal positions ����ڵ�λ�� sevec3�������ͣ���vectorsimdͷ�ļ��ж���

	//==== input mesh topology data for collision computation ����������������������ײ����

	const Int4* m_edges;				// indices of the two adjacent vertices and the two opposite vertices of edges �����������ڶ������Զ���ıߵ���������ά������ÿһά��������
	const Int4* m_faces;				// indices of the three adjacent vertices of faces, the fourth is useless �������������ڶ��㣬���ĸ�λ����Ч
	
	const SeCsr<int>* m_neighbours;	// the adjacent information of vertices stored in a csr format ��csr��ʽ�洢������ڽӾ���

public:
	SeSchwarzPreconditioner() {};
	//==== call before time integration once a frame ʱ�����ǰ����һ��  ÿһ֡��ʱ�����
	void AllocatePrecoditioner(int numVerts, int numEdges, int numFaces);

	//==== call before PCG iteration loop PCG����ѭ��ǰ����
	void PreparePreconditioner(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges);

	//void PreparePreconditioner(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges,
	//	const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);


	//==== call during PCG iterations  PCG�����ڼ����
	void Preconditioning(SeVec3fSimd* z, const SeVec3fSimd* residual, int dim); //r��ʼȡb
	std::vector<int>                    m_MapperSortedGetOriginal;
	int m_numVerts = 0;

private:
	// �������� + ��������; Ҫ�����ױ����ͺ��������͡�����
	 //����������߸���&�����θ���
	int m_numEdges = 0;
	int m_numFaces = 0;

	int m_frameIndex = 0; //֡����

	int m_totalSz = 0; ////ÿһ����صĸ���*size(32)*1.5 ��Ȼ���ۼӵĽ��
	int m_numLevel = 0; //�༶level�ļ���

	int m_totalNumberClusters; //total��

	
	SeAabb<SeVec3fSimd>					m_aabb;


	SeArray2D<SeMatrix3f>				m_hessianMapped;  //����hessianӳ��
	std::vector<int>					m_mappedNeighborsNum; //Ī�ٿռ���ÿ��������ھӸ���(�����Լ�)
	std::vector<int>					m_mappedNeighborsNumRemain; //���������ʹ�õ�remain����
	SeArray2D<int>						m_mappedNeighbors;  // ��������󶥵���ھ������б�
	SeArray2D<int>						m_mappedNeighborsRemain;

	SeArray2D<int>						m_CoarseSpaceTables; //ϵ���ռ���
	std::vector<int>                    m_prefixOrignal; //��ʼǰ׺
	std::vector<unsigned int>           m_fineConnectMask; 
	std::vector<unsigned int>           m_nextConnectMsk;
	std::vector<unsigned int>           m_nextPrefix; //

	std::vector<SeVec3fSimd>			m_mappedPos; //Ī����ӳ��������

	std::vector<Int4>					m_coarseTables;
	//goingnext:�����б�,��ʾ��ǰ��������һ����������ӳ�� m_goingNext[0]=k��ʾ��ǰlevel����Ϊ0�ĵ�����һ��level������Ϊk;���Ⱥܳ�������������level�����ӳ�䣺
	//ע����һ��level��Ų��ٴ�0��ʼ����Ҫ����һ�����ж���ĸ�����Ϊƫ����
	std::vector<int>					m_goingNext; 
	std::vector<int>					m_denseLevel;
	std::vector<Int2>					m_levelSize;  // .x = current level size ��ǰlevel��ʵ�ʶ������(����������Ķ���-����һ����ǿ�д�������)   .y = (prefixed) current level begin index (ǰ׺)��ǰ����ʼ��������ǰlevel֮ǰ����level�Ķ���ĸ�����������ǰ(�����ǰ������ⶥ���)

	std::vector<SeVec3fSimd>			m_mappedR;
	std::vector<SeVec3fSimd>			m_mappedZ;
	//std::vector<int>                    m_MapperSortedGetOriginal;          // sorted by morton 
	std::vector<int>                    m_mapperOriginalGetSorted;   
	std::vector<SeMorton64>             m_mortonCode;  //Ī�ٱ���

	SeArray2D<SeMatrix3f>				m_hessian32;   //�洢����hessianϵͳ���� = ��������ײ(��ʼ����)+������ײ����hessian
	std::vector<SeMatrix3f>             m_additionalHessian32; //��3*3���󹹳ɵ�����
	std::vector<float>                  m_invSymR; //ϵͳ�����

	
	int m_stencilNum;
	int m_maxStencilNum;
	std::vector<Stencil>				m_stencils;  
	std::vector<Int5>					m_stencilIndexMapped;

private:

	void ComputeLevelNums(int bankSize); //���㼶�����

	void DoAlllocation(); //��λ��
	 
	void ComputeTotalAABB(); //��AABB����

	void ComputeAABB(); //��AABB

	void SpaceSort(); //�ռ�����(Ī��������)

	void FillSortingData();//����������� 

	void DoingSort(); //����

	void ComputeInverseMapper(); //����

	void MapHessianTable();//���׾���ӳ��

	void PrepareCollisionStencils(const EfSet* efSets, const EeSet* eeSets, const VfSet* vfSets, unsigned int* efCounts, unsigned int* eeCounts, unsigned int* vfCounts);

	void MapCollisionStencilIndices(); //��ײ��һЩ����

	void ReorderRealtime(); //ʵʱ��������



	void BuildConnectMaskL0(); //����ӳ��L0

	void BuildCollisionConnection(unsigned int* pConnect, const int* pCoarseTable);

	void PreparePrefixSumL0();//

	void BuildLevel1(); //L1

	void BuildConnectMaskLx(int level); //����LX

	void NextLevelCluster(int level); //��һ���Ĵ�

	void PrefixSumLx(int level); //sum

	void ComputeNextLevel(int level);//������һ��

	void TotalNodes();//����ڵ���

	void AggregationKernel();//�����ں�

	void AdditionalSchwarzHessian2(SeMatrix3f hessian, std::vector<SeMatrix3f>& pAdditionalHessian, SeArray2D<SeMatrix3f>& pDenseHessian, int v1, int v2, const std::vector<int>& pGoingNext, int nLevel, int vNum);

	void PrepareCollisionHessian(); //��ײ׼��

	void PrepareHessian(const SeMatrix3f* diagonal, const SeMatrix3f* csrOffDiagonals, const int* csrRanges); //���׾���׼��
	
	void LDLtInverse512();//LDLT 

	void BuildResidualHierarchy(const SeVec3fSimd* m_cgResidual);

	void SchwarzLocalXSym();

	void CollectFinalZ(SeVec3fSimd* m_cgZ); //����ռ�
};

SE_NAMESPACE_END