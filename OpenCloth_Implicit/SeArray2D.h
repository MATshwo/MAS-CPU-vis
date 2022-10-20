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

#include "SePreDefine.h"

SE_NAMESPACE_BEGIN

/*************************************************************************
****************************    SeArray2D    *****************************
*************************************************************************/

//!	@brief	二维数组容器——这个二维数组包含三个属性，列、行、还有一个(按列)存储二维数组值的列向量
// 目前来看，定义该类的目的会有助于稀疏结构的存储
//这个类的作用是将一个二维矩阵按行/列的形式铺平成一维向量的形式进行存储，里边的type一般是int/float

template<typename Type> class SeArray2D
{

public:

	//!	@brief	Default constructor.
	SeArray2D() : m_Rows(0), m_Columns(0) {}

	//!	@brief	Construct and resize.默认给向量m_value赋值为0
	explicit SeArray2D(size_t _Rows, size_t _Columns) : m_Rows(_Rows), m_Columns(_Columns), m_Values(_Rows * _Columns) {} 

	//!	@brief	Construct, resize and memory set.
	explicit SeArray2D(size_t _Rows, size_t _Columns, Type _Value) : m_Rows(_Rows), m_Columns(_Columns), m_Values(_Rows * _Columns, _Value) {}

public:

	//! @brief  Data will be reserved in Column Major format 列向量存储
	void Resize(size_t _Rows, size_t _Columns)
	{
		m_Values.resize(_Rows * _Columns);

		m_Columns = _Columns;

		m_Rows = _Rows;
	}

	//! @brief  release spare memory
	void ShrinkToFit()
	{
		m_Values.shrink_to_fit(); //释放空闲内存，让size&capacity保持一致
	}

	//!	@brief	Fill data. 数据填充，这里每个元素都用同一个值填充
	void Memset(Type _Value, size_t _Begin = 0, size_t _End = SIZE_MAX)
	{
		_End = _End < m_Values.size() ? _End : m_Values.size(); //指定size进行填充

		for (size_t i = _Begin; i < _End; ++i)
		{
			m_Values[i] = _Value;
		}
	}

	//!	@brief	Exchange context with right. 
	void Swap(SeArray2D & _Right)
	{
		m_Values.swap(_Right.m_Values); //向量交换：size和capacity和value
		 
		SE_SWAP(m_Columns, _Right.m_Columns); //对应结构体内的行列信息也交换更新

		SE_SWAP(m_Rows, _Right.m_Rows);
	}

	void Clear()
	{
		m_Values.clear();

		m_Values.shrink_to_fit(); //向量元素清除，并释放内存，行列信息也清除=0

		m_Rows = m_Columns = 0;
	}

public:

	const Type * operator[](size_t i) const { SE_ASSERT(i < m_Rows);  return &m_Values[m_Columns * i]; } //[]运算符重载，返回第i行元素的地址  A[i][j]对应i行j列的元素，也就是向量第i*m_columns+j个元素

	Type * operator[](size_t i) { SE_ASSERT(i < m_Rows);  return &m_Values[m_Columns * i]; }

	const Type * Ptr() const { return m_Values.data(); } //返回第一个元素值(&地址)，不支持修改

	bool IsEmpty() const { return m_Values.empty(); } //判断是否为空

	size_t Size() const { return m_Values.size(); }//返回向量s长度

	size_t Capacity() const { return m_Values.capacity(); } //返回已分配内存大小，也就是该结构体的大小容量=向量存储的大小和容量

	size_t Columns() const { return m_Columns; }

	size_t Rows() const { return m_Rows; }

	Type * Ptr() { return m_Values.data(); } //返回当前向量的第一个元素地址，且支持修改:实质就是返回m_values这个向量

private:

	std::vector<Type> m_Values; //一维向量，元素类型为type类型，列向量的形式

	size_t m_Rows, m_Columns;   //无符号整型的别名 
};

SE_NAMESPACE_END