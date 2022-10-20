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

#include "SeMathSimd.h"

#include "SeAabb.h"

#include <limits>

SE_NAMESPACE_BEGIN

/*************************************************************************
******************************    SeAabb    ******************************
*************************************************************************/
 
/**
 *	@brief		Axis-align Bounding Box. //simd？好像是为伪四维向量特殊设置的存在，vector\math代码块中均有该部分存在
 */

//SeAabb.h定义了模板结构体，指定输入一个任意格式的数据类型，并要求该数据类型里元素的取值是浮点型；
//而此处则明确指定结构体的输入类别是SimdVector3<float>; ——即长度为4，实质有效数值为3的存储顶点数据的特殊向量的类型。

template<> struct SeAabb<Simd3f>
{

public:

	//!	@brief	Default constructor.
	SeAabb() : Lower(FLT_MAX), Upper(-FLT_MAX) {} //初始包围盒无穷大

	//!	@brief	Constructed by a given point.
	SeAabb(const Simd3f & Point) : Lower(Point), Upper(Point) {}

	//!	@brief	Constructed by a paired points explicitly.
	explicit SeAabb(const Simd3f & P1, const Simd3f & P2) : Lower(Math::Min(P1, P2)), Upper(Math::Max(P1, P2)) {}

public:

	SeAabb operator+(const SeAabb & Other) const
	{
		return SeAabb(Math::Min(Lower, Other.Lower), Math::Max(Upper, Other.Upper));
	}

	SeAabb operator+(const Simd3f & Point) const
	{
		return SeAabb(Math::Min(Lower, Point), Math::Max(Upper, Point));
	}

	void operator += (const SeAabb & Other)
	{
		Lower = Math::Min(Lower, Other.Lower);		Upper = Math::Max(Upper, Other.Upper);
	}

	void operator += (const Simd3f & Point)
	{
		Lower = Math::Min(Lower, Point);		Upper = Math::Max(Upper, Point);
	}



	void Enlarge(float Amount)
	{
		Lower -= Amount;		Upper += Amount;
	}

	Simd3f Center() const
	{
		return Math::Average(Lower, Upper);
	}

	Simd3f Extent() const
	{
		return Upper - Lower;
	}

public:

	Simd3f Lower, Upper;
};


//是否包含：顶点p是否包含在aabb这个包围盒内
template<>
static SE_INLINE bool IsContain(const SeAabb<Simd3f>& aabb, const Simd3f& P) 
{ 
	return !((aabb.Upper.x < P.x) || (aabb.Upper.y < P.y) || (aabb.Upper.z < P.z) ||
			 (aabb.Lower.x > P.x) || (aabb.Lower.y > P.y) || (aabb.Lower.z > P.z));
}

//此处的代码是符合直观的，这里也说明了aabb.h这里的代码确实存在问题
//它想实现的功能就是检测顶点P是否包含在扩充radium的aabb内部（条件放宽了）
template<> 
static SE_INLINE bool IsContain(const SeAabb<Simd3f>& aabb, const Simd3f& P, float radium)
{ 
	return !((aabb.Upper.x < P.x - radium) || (aabb.Upper.y < P.y - radium) || (aabb.Upper.z < P.z - radium) ||
		(aabb.Lower.x > P.x + radium) || (aabb.Lower.y > P.y + radium) || (aabb.Lower.z > P.z + radium));
}

//是否重叠，返回1表示重叠
template<>
static SE_INLINE bool IsOverlap(const SeAabb<Simd3f>& aabb, const SeAabb<Simd3f>& P)
{
	return !((aabb.Upper.x < P.Lower.x) || (aabb.Upper.y < P.Lower.y) || (aabb.Upper.z < P.Lower.z) ||
		(aabb.Lower.x > P.Upper.x) || (aabb.Lower.y > P.Upper.y) || (aabb.Lower.z > P.Upper.z));
}

//是否相交
template<>
static SE_INLINE bool IsIntersect(const SeAabb<Simd3f>& aabb, const Simd3f& Pa, const Simd3f& Pb) 
{
	Simd3f Dir = Simd3f(Pb - Pa);

	if (Dir.x == 0.0f)		Dir.x = 1e-6f;
	if (Dir.y == 0.0f)		Dir.y = 1e-6f;
	if (Dir.z == 0.0f)		Dir.z = 1e-6f;

	Simd3f invDir = 1.0f / Dir;
	Simd3f Left = Simd3f(aabb.Lower - Pa) * invDir;
	Simd3f Right = Simd3f(aabb.Upper - Pa) * invDir;

	Float2 T1 = Float2(Math::Min(Left.x, Right.x), Math::Max(Left.x, Right.x));
	Float2 T2 = Float2(Math::Min(Left.y, Right.y), Math::Max(Left.y, Right.y));
	Float2 T3 = Float2(Math::Min(Left.z, Right.z), Math::Max(Left.z, Right.z));

	float range0 = Math::Max(T1.x, T2.x, T3.x, 0.0f);
	float range1 = Math::Min(T1.y, T2.y, T3.y, 1.0f);

	return range0 < range1;
}

SE_NAMESPACE_END
