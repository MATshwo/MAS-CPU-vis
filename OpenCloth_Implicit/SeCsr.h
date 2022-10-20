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

//压缩稀疏数据——将已有的稀疏矩阵信息提取
template<typename Type> class SeCompressSparseData
{
public:

	SeCompressSparseData() {}
	SeCompressSparseData(const std::vector<int>& starts, const std::vector<int>& idxs, const std::vector<Type>& values) //三个指标实施稀疏数据存储:start(每行首个非零元的位置),idx(非零索引) and value(存储非零值)
		:
		m_starts(starts), //整数向量：每一行的开始索引
		m_idxs(idxs),  //整数向量：对应列的索引
		m_values(values)  //指定类型的向量：存储非零元素
	{}

	virtual ~SeCompressSparseData() {}

	//考虑代码在主文件中的位置，这里的初始化idx执行的具体功能是什么？只是初始化start和idx，并未对value赋值
	//如果要稀疏数据压缩的功能:输入的数组需要有明确的含义：1）array2D格式，向量的每个元素是对应行非零元素的列索引;2）SeArray2D格式，arr[i][j]存储的数据必须满足，arr[i][0]存储该行元素的个数，剩余存储该行非零元素的列索引。
	//满足上述条件输入，对应下列函数才能执行到压缩数据的功能:提取行数组+列索引数组


	void InitIdxs(const std::vector<std::vector<int>>& array2D) //二维矩阵：元素为向量的向量
	{
		Clear();
		//输入：[[3,2,1],[1,1,1]],size=2,m_start=[0,3,3+3],m_idsx=[3,2,1,1,1,1]
		int size = array2D.size();

		m_starts.resize(size + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < size; ++i)
		{
			int num = array2D[i].size(); //结合CSR的存储格式，num表示非零元素的个数，但此处使用了size(包含0元素和非零元素)，难以理解！

			m_starts[i + 1] = m_starts[i] + num;
		}

		m_idxs.resize(m_starts.back()); //总元素的个数，实验了下就算没有初始化m_id的值运行之后也会初始化为0，所以没有具体id值指向时默认=0；

		OMP_PARALLEL_FOR //开启多线程
			for (int i = 0; i < size; ++i)
			{
				int num = array2D[i].size(); //m_idxs.data()表示数组首元素的地址;+m[i]表示原始地址偏移m_starts[i]个type类型长度 (例如，m[0]=3占据3*4(=sizeof(int))=12字节，就在首地址基础上+12)
				std::memcpy(m_idxs.data() + m_starts[i], array2D[i].data(), sizeof(int) * num); //赋值操作：将array2D中的数据按索引赋值,
			}
	}
	void InitIdxs(const SeArray2D<int>& array2D, bool isRowMajor)
	{
		//输入数据类型是array2D
		Clear();
		//输入arr.value = [2,1,3,2,1,1,2,7,1] r*c=3*3，row=3,m_start=[0,2,4,6],m_idx=[]
		int DIM = isRowMajor ? array2D.Rows() : array2D.Columns();

		m_starts.resize(DIM + 1); Utility::MemsetZero(m_starts);

		for (int i = 0; i < DIM; ++i)
		{
			int num = isRowMajor ? array2D[i][0] : array2D[0][i]; //[]取地址+隐式格式转换，对应取i行0列的值

			//此处也有上方不对应，这里num并非表示每一行的长度，只是记录了每行首元素的值
			//结合address[k - 1] = isRowMajor ? array2D[i][k] : array2D[k][i]; 这一代码大致可以推断出该数组类型每行的首元素就是该行元素的个数，这样才能与上文对应
			//推断不一定正确,但这样首先可以将上下两种压缩方法对应，表达同一个意思
			//其次需要考虑如何体现出压缩数据的思想：目前不考虑特殊输入格式前提下，该代码并未有效区分零元素&非零元素，只能寄希望于输入的数组/2D是已经处理的形式，具体形式需要结合上下文
			m_starts[i + 1] = m_starts[i] + num; //同上，记录元素个数
		}
		m_idxs.resize(m_starts.back());

		OMP_PARALLEL_FOR
			for (int i = 0; i < DIM; ++i)
			{
				int* address = m_idxs.data() + m_starts[i]; //int *p = a[0]+3;表示p的地址指向a首元素偏移三个sizes,这里sizes=sizeof(a[0]);但注意：
				// 如果int *p = (int*)(&a+1);要注意这里偏移的sizes = sizeof(a)而不是一个元素的长度； sizeof(p)=sizeof(p的指针类型)
				int num = Size(i); //num对应第i行有多少元素

				for (int k = 1; k <= num; ++k)
				{
					address[k - 1] = isRowMajor ? array2D[i][k] : array2D[k][i]; //第i行k列的元素
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
		return m_starts.back(); //总共有多少非零元素
	}

	int Size(int id) const
	{
		return -m_starts[id] + m_starts[id + 1]; //相邻行/列之间有几个元素 ：按照稀疏矩阵压缩算法，这里是指id-1行有多少非0元素
	}

	int Start(int id) const
	{
		return m_starts[id]; //前id-1行有多少非零元素
	}

	const int* StartPtr(int id) const
	{
		return &(m_starts[id]); //取地址
	}

	const int Idx(int id) const
	{
		return m_idxs[id]; //第id个元素的里下标
	}

	const int* IdxPtr(int id) const
	{
		return m_idxs.data() + m_starts[id]; //返回id行m_idx首元素的地址，行逻辑访问该行非零元素的列索引
	}

	const Type* ValuePtr(int id) const
	{
		return m_values.data() + m_starts[id]; //id行m_value首元素地址，行逻辑访问非零元素值
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

	int Rows() const { return SeCompressSparseData<Type>::m_starts.size() - 1; } //行数=start-1

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
		return this; //返回当前对象的地址
	}
};

SE_NAMESPACE_END

