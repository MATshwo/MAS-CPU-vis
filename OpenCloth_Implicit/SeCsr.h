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

#include "SeOmp.h"
#include "SeArray2D.h"
#include "SeUtility.h"

SE_NAMESPACE_BEGIN

//ѹ��ϡ�����ݡ��������е�ϡ�������Ϣ��ȡ
template<typename Type> class SeCompressSparseData
{
public:

	SeCompressSparseData() {}
	SeCompressSparseData(const std::vector<int>& starts, const std::vector<int>& idxs, const std::vector<Type>& values) //����ָ��ʵʩϡ�����ݴ洢:start(ÿ���׸�����Ԫ��λ��),idx(��������) and value(�洢����ֵ)
		:
		m_starts(starts), //����������ÿһ�еĿ�ʼ����
		m_idxs(idxs),  //������������Ӧ�е�����
		m_values(values)  //ָ�����͵��������洢����Ԫ��
	{}

	virtual ~SeCompressSparseData() {}

	//���Ǵ��������ļ��е�λ�ã�����ĳ�ʼ��idxִ�еľ��幦����ʲô��ֻ�ǳ�ʼ��start��idx����δ��value��ֵ
	//���Ҫϡ������ѹ���Ĺ���:�����������Ҫ����ȷ�ĺ��壺1��array2D��ʽ��������ÿ��Ԫ���Ƕ�Ӧ�з���Ԫ�ص�������;2��SeArray2D��ʽ��arr[i][j]�洢�����ݱ������㣬arr[i][0]�洢����Ԫ�صĸ�����ʣ��洢���з���Ԫ�ص���������
	//���������������룬��Ӧ���к�������ִ�е�ѹ�����ݵĹ���:��ȡ������+����������


	void InitIdxs(const std::vector<std::vector<int>>& array2D) //��ά����Ԫ��Ϊ����������
	{
		Clear();
		//���룺[[3,2,1],[1,1,1]],size=2,m_start=[0,3,3+3],m_idsx=[3,2,1,1,1,1]
		int size = array2D.size();

		m_starts.resize(size + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < size; ++i)
		{
			int num = array2D[i].size(); //���CSR�Ĵ洢��ʽ��num��ʾ����Ԫ�صĸ��������˴�ʹ����size(����0Ԫ�غͷ���Ԫ��)��������⣡

			m_starts[i + 1] = m_starts[i] + num;
		}

		m_idxs.resize(m_starts.back()); //��Ԫ�صĸ�����ʵ�����¾���û�г�ʼ��m_id��ֵ����֮��Ҳ���ʼ��Ϊ0������û�о���idֵָ��ʱĬ��=0��

		OMP_PARALLEL_FOR //�������߳�
			for (int i = 0; i < size; ++i)
			{
				int num = array2D[i].size(); //m_idxs.data()��ʾ������Ԫ�صĵ�ַ;+m[i]��ʾԭʼ��ַƫ��m_starts[i]��type���ͳ��� (���磬m[0]=3ռ��3*4(=sizeof(int))=12�ֽڣ������׵�ַ������+12)
				std::memcpy(m_idxs.data() + m_starts[i], array2D[i].data(), sizeof(int) * num); //��ֵ��������array2D�е����ݰ�������ֵ,
			}
	}
	void InitIdxs(const SeArray2D<int>& array2D, bool isRowMajor)
	{
		//��������������array2D
		Clear();
		//����arr.value = [2,1,3,2,1,1,2,7,1] r*c=3*3��row=3,m_start=[0,2,4,6],m_idx=[]
		int DIM = isRowMajor ? array2D.Rows() : array2D.Columns();

		m_starts.resize(DIM + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < DIM; ++i)
		{
			int num = isRowMajor ? array2D[i][0] : array2D[0][i]; //[]ȡ��ַ+��ʽ��ʽת������Ӧȡi��0�е�ֵ

			//�˴�Ҳ���Ϸ�����Ӧ������num���Ǳ�ʾÿһ�еĳ��ȣ�ֻ�Ǽ�¼��ÿ����Ԫ�ص�ֵ
			//���address[k - 1] = isRowMajor ? array2D[i][k] : array2D[k][i]; ��һ������¿����ƶϳ�����������ÿ�е���Ԫ�ؾ��Ǹ���Ԫ�صĸ������������������Ķ�Ӧ
			//�ƶϲ�һ����ȷ,���������ȿ��Խ���������ѹ��������Ӧ�����ͬһ����˼
			//�����Ҫ����������ֳ�ѹ�����ݵ�˼�룺Ŀǰ���������������ʽǰ���£��ô��벢δ��Ч������Ԫ��&����Ԫ�أ�ֻ�ܼ�ϣ�������������/2D���Ѿ��������ʽ��������ʽ��Ҫ���������
			m_starts[i + 1] = m_starts[i] + num; //ͬ�ϣ���¼Ԫ�ظ���
		}
		m_idxs.resize(m_starts.back());

		OMP_PARALLEL_FOR
			for (int i = 0; i < DIM; ++i)
			{
				int* address = m_idxs.data() + m_starts[i]; //int *p = a[0]+3;��ʾp�ĵ�ַָ��a��Ԫ��ƫ������sizes,����sizes=sizeof(a[0]);��ע�⣺
				// ���int *p = (int*)(&a+1);Ҫע������ƫ�Ƶ�sizes = sizeof(a)������һ��Ԫ�صĳ��ȣ� sizeof(p)=sizeof(p��ָ������)
				int num = Size(i); //num��Ӧ��i���ж���Ԫ��

				for (int k = 1; k <= num; ++k)
				{
					address[k - 1] = isRowMajor ? array2D[i][k] : array2D[k][i]; //��i��k�е�Ԫ��
				}
			}
	}
	void Clear()
	{
		Utility::ClearAndShrink(m_starts);
		Utility::ClearAndShrink(m_idxs);
		Utility::ClearAndShrink(m_values);
	}

	int Size() const
	{
		return m_starts.back(); //�ܹ��ж��ٷ���Ԫ��
	}

	int Size(int id) const
	{
		return -m_starts[id] + m_starts[id + 1]; //������/��֮���м���Ԫ�� ������ϡ�����ѹ���㷨��������ָid-1���ж��ٷ�0Ԫ��
	}

	int Start(int id) const
	{
		return m_starts[id]; //ǰid-1���ж��ٷ���Ԫ��
	}

	const int* StartPtr(int id) const
	{
		return &(m_starts[id]); //ȡ��ַ
	}

	const int Idx(int id) const
	{
		return m_idxs[id]; //��id��Ԫ�ص����±�
	}

	const int* IdxPtr(int id) const
	{
		return m_idxs.data() + m_starts[id]; //����id��m_idx��Ԫ�صĵ�ַ�����߼����ʸ��з���Ԫ�ص�������
	}

	const Type* ValuePtr(int id) const
	{
		return m_values.data() + m_starts[id]; //id��m_value��Ԫ�ص�ַ�����߼����ʷ���Ԫ��ֵ
	}

	virtual const SeCompressSparseData* Ptr() const
	{
		return this;
	}

protected:

	std::vector<int>	m_starts;
	std::vector<int>	m_idxs;
	std::vector<Type>	m_values;
};

template<typename Type> class SeCsr : public SeCompressSparseData<Type>
{
public:
	SeCsr() {}
	SeCsr(const std::vector<int>& starts, const std::vector<int>& idxs, const std::vector<Type>& values) :SeCompressSparseData<Type>(starts, idxs, values) {}

	int Rows() const { return SeCompressSparseData<Type>::m_starts.size() - 1; } //����=start-1

	virtual const SeCsr* Ptr() const override
	{
		return this;
	}
};


template<typename Type> class SeCsc : public SeCompressSparseData<Type>
{
public:

	int Columns() const { return SeCompressSparseData<Type>::m_starts.size() - 1; }

	virtual const SeCsc* Ptr() const override
	{
		return this; //���ص�ǰ����ĵ�ַ
	}
};

SE_NAMESPACE_END

