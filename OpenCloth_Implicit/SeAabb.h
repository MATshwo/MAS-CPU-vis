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

#include "SeMath.h"
#include <cfloat>

SE_NAMESPACE_BEGIN

/*************************************************************************
******************************    SeAabb    ******************************
*************************************************************************/
 
/**
 *	@brief		Axis-align Bounding Box. 构造轴对齐包围盒
 */
//注意此结构体的含义：此结构体对应一个名称为包围盒的数据类型，接受一个向量或多个向量作为输入，具有特定的上下界属性，能够通过.lower访问其对应包围盒的上下界。aabb是这个包围盒的名称，它具有上下界以及加减运算等一系列的属性。
template<typename VecType> struct SeAabb
{
	//这里的包围盒由于可以对向量比较大小，所以可以返回一个具体包围盒的详细参数：长宽高...
	// 仅支持浮点型数据？不理解：再看明白了，是指输入向量对应元素的格式必须是浮点型！
	static_assert(std::is_floating_point<typename VecType::value_type>::value, "Only support floating type currently.");

public:

	//!	@brief	Default constructor. 没有数据输入时使用C++中的常量进行赋值
	 SeAabb() : Lower(FLT_MAX), Upper(-FLT_MAX) {} //下界和上界（正负无穷大）

	//!	@brief	Constructed by a given point. 传入一个参数的情况
	 SeAabb(const VecType & Point) : Lower(Point), Upper(Point) {}

	//!	@brief	Constructed by a paired points explicitly. 给定一对数据进行赋值
	 explicit SeAabb(const VecType & P1, const VecType & P2) : Lower(Math::Min(P1, P2)), Upper(Math::Max(P1, P2)) {}

	//!	@brief	Constructed by a paired points explicitly. 三个值进行赋值
	//牛！！果然大佬都是自己编写数学库的：SeMath头文件中的max/min函数，可以对标量向量比较大小，如果是向量，返回每个维度取最大/小值后得到的新的向量。
	 explicit SeAabb(const VecType & P1, const VecType & P2, const VecType & P3) : Lower(Math::Min(P1, P2, P3)), Upper(Math::Max(P1, P2, P3)) {}

public:

	//运算符重载: 根据传入参数不同分开定义；
	 SeAabb<VecType> operator+(const SeAabb<VecType> & Other) const
	{
		return SeAabb<VecType>(Math::Min(Lower, Other.Lower), Math::Max(Upper, Other.Upper));
	}

	 SeAabb<VecType> operator+(const VecType & Point) const
	{
		return SeAabb<VecType>(Math::Min(Lower, Point), Math::Max(Upper, Point));
	}

	 void operator+=(const SeAabb<VecType> & Other)
	{
		Lower = Math::Min(Lower, Other.Lower);		Upper = Math::Max(Upper, Other.Upper);
	}

	 void operator+=(const VecType & Point)
	{
		Lower = Math::Min(Lower, Point);		Upper = Math::Max(Upper, Point);
	}

	 

	 void Enlarge(float Amount) //拓宽范围，上下各扩展amount
	{
		Lower -= Amount;		Upper += Amount;
	}

	 VecType Center() const //求中心位置
	{
		return Math::Average(Lower, Upper);
	}

	 VecType Extent() const //范围
	{
		return Upper - Lower;
	}

public:
	// 搞不懂为什么是在最后声明变量，是有什么规定？
	VecType Lower, Upper;
};


//以下包含三类函数：功能为了检测包含、相交、重合的情况
/*************************************************************************
***********************    IsInside<AABB-AABB>   ************************ 
 *************************************************************************/
//检测是否有包含

template<typename VecT> 
static SE_INLINE bool IsContain(const SeAabb<VecT>& aabb, const VecT& P) 
{ 
	//是否包含：检测P是否在aabb对应包围盒的内部
	const int componentCount = VecT::Elements;
	for (int i = 0; i < componentCount; i++)
	{
		if (P[i]<aabb.Lower[i] || P[i]>aabb.Upper[i])
		{
			//不在aabb对应包围盒的内部
			return false;
		}
	}
	return true; 

}
template<typename VecT> 
static SE_INLINE bool IsContain(const SeAabb<VecT>& aabb, const VecT& P,float radium) // static __inline 静态内联函数：只在当前文件内有效的内联函数
{ 
	//扩大后的aabb是否还在aabb内部？？？这里好像没有用到P点,总感觉怪怪的，这不应该肯定是false,可能是负数？
	//感觉这个地方有点问题，这个应该是为了检测P是否在aabb扩充的包围盒内部，所以return 应该是enlargedAabb和P？
	auto enlargedAabb = aabb; 
	enlargedAabb.Upper += radium; //扩充aabb包围盒
	enlargedAabb.Lower -= radium;
	return IsContain(aabb, enlargedAabb); 
}

template<> 
static SE_INLINE bool IsContain(const SeAabb<Float2>& aabb, const Float2& P) //指定二维浮点型
{ 
	return !((aabb.Upper.x < P.x) || (aabb.Upper.y < P.y) ||
			 (aabb.Lower.x > P.x) || (aabb.Lower.y > P.y));
}

template<> 
static SE_INLINE bool IsContain(const SeAabb<Float3>& aabb, const Float3& P) 
{ 
	return !((aabb.Upper.x < P.x) || (aabb.Upper.y < P.y) || (aabb.Upper.z < P.z) ||
			 (aabb.Lower.x > P.x) || (aabb.Lower.y > P.y) || (aabb.Lower.z > P.z));
}

/*************************************************************************
***********************    IsInside<AABB-Line>    ************************
*************************************************************************/
//检测包围盒和线条的关系：是否有相交的情况
template<typename VecT>
static  bool IsIntersect(const SeAabb<VecT>& aabb, const VecT& Pa, const VecT& Pb, const VecT& iv);

template<typename VecT>
static  bool IsIntersect(const SeAabb<VecT>& aabb, const VecT& Pa, const VecT& Pb);


template<>
static SE_INLINE bool IsIntersect(const SeAabb<Float3>& aabb, const Float3& Pa, const Float3& Pb) 
{
	Float3 Dir = Float3(Pb - Pa); //检测pa&pb的连线

	if (Dir.x == 0.0f)		Dir.x = 1e-6f;
	if (Dir.y == 0.0f)		Dir.y = 1e-6f;
	if (Dir.z == 0.0f)		Dir.z = 1e-6f;

	Float3 invDir = 1.0f / Dir;
	Float3 Left = Float3(aabb.Lower - Pa) * invDir; //以二维为例：这里实质使用上下界到左顶点长度除以线段长度，根据相交/不相交的情况可以得出结论
	Float3 Right = Float3(aabb.Upper - Pa) * invDir;

	Float2 T1 = Float2(Math::Min(Left.x, Right.x), Math::Max(Left.x, Right.x)); //x方向线条在包围盒左侧还是右侧，返回二维向量指
	Float2 T2 = Float2(Math::Min(Left.y, Right.y), Math::Max(Left.y, Right.y)); //y...
	Float2 T3 = Float2(Math::Min(Left.z, Right.z), Math::Max(Left.z, Right.z)); //z...

	float range0 = Math::Max(T1.x, T2.x, T3.x, 0.0f); //
	float range1 = Math::Min(T1.y, T2.y, T3.y, 1.0f); //

	return range0 < range1; //true表示相交，画了图大致知道是什么意思
}

//检测是否重叠 ，直接false可还行
template<typename VecT>
static SE_INLINE bool IsOverlap(const SeAabb<VecT>& aabb, const SeAabb<VecT>& P) { return false; }


SE_NAMESPACE_END
