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

namespace Utility
{
	//感觉在这个空间内定义了一系列给data向量赋值的方法：一个向量赋值、常量赋值、不同数据类型的向量赋值、。。。
	// size_t = unsigned int 
	//有的类似赋值会出现两种定义：主要跟参数类型不同，一个是以type为元素的向量，一个是type类型的变量
	template <typename Type>
	size_t ByteSize(const std::vector<Type> & data)
	{
		return sizeof(Type) * data.size(); //返回data向量总长度：=向量长度*每个元素的长度
	}

	//使用src变量赋值
	template <typename Type>
	void Memcpy(Type * dstData, const Type * srcData, size_t size)
	{
		std::memcpy(dstData, srcData, size * sizeof(Type)); //内存拷贝函数，第三个参数表示复制的长度
	}

	//使用另一种type的向量赋值，需要进行类型转换，对长的那个进行裁剪
	template <typename Type1, typename Type2>
	void Memcpy(std::vector<Type1> & data, const std::vector<Type2> & srcData, size_t size)
	{
		/*if(ByteSize(data) < ByteSize(srcData))
		{
			SE_WARNING_LOG("DesSize is less than Src Size!");
		}*/
		std::memcpy(data.data(), srcData.data(), SE_MIN(SE_MIN(ByteSize(data), ByteSize(srcData)), size * sizeof(Type1)));
	}

	template <typename Type>
	void Memset(std::vector<Type>& data, const Type & value)
	{
		std::fill(data.begin(), data.end(), value); //用value的值填充data的每个元素
	}

	template <typename Type>
	void Memset(Type* data, const Type& value, size_t size)
	{
		std::fill(data, data + size, value);
	}

	//全0赋值
	template <typename Type>
	void MemsetZero(std::vector<Type> & data)
	{
		std::memset(data.data(), 0, ByteSize(data));
	}

	template <typename Type>
	void MemsetZero(Type * data, size_t size)
	{
		std::memset(data, 0, sizeof(Type) * size);
	}

	//-1赋值
	template <typename Type>
	void MemsetMinusOne(std::vector<Type> & data)
	{
		std::memset(data.data(), -1, ByteSize(data));
	}

	template <typename Type>
	void MemsetMinusOne(Type * data, size_t size)
	{
		std::memset(data, -1, sizeof(Type) * size);
	}

	//清除操作：减少容器的容量以适应其大小并销毁超出容量的所有元素。
	template <typename Type>
	void ClearAndShrink(std::vector<Type> & data)
	{
		data.clear(); data.shrink_to_fit();
	}

	//resize操作
	template <typename Type>
	void ResizeAndShrink(std::vector<Type> & data, size_t dim)
	{
		data.resize(dim); data.shrink_to_fit();
	}

	//copy操作
	template <typename Type>
	void CopyAndShrink(std::vector<Type> & desData, const std::vector<Type> & srcData)
	{
		desData = srcData; desData.shrink_to_fit();
	}
}
SE_NAMESPACE_END
