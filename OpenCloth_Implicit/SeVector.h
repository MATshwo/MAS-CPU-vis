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

//自定义的vector结构——这里定义了多个向量结构体并对其类型进行命名，本框架中所有变量类型基本都参照本部分；
#include "SePreDefine.h"

#include <cmath>

SE_NAMESPACE_BEGIN

#define SE_VEC2_ALIGN(Type)		SE_ALIGN(SE_MIN(sizeof(Type) * 2, 16))		//!	Maximum align bits -> 64.
#define SE_VEC4_ALIGN(Type)		SE_ALIGN(SE_MAX(sizeof(Type) * 4, 16))		//!	Maximum align bits -> 64.

/*************************************************************************
*****************************    SDVector    *****************************
*************************************************************************/
//枚举类型
enum VecDescription : int { SIMD_VEC3 = -1, VEC2_XY = 2, VEC3_XYZ = 3, VEC4_XYZW = 4 };

template<typename Type, int N> struct SeVector;
template<typename Type> using SeVector2 = SeVector<Type, VEC2_XY>; //2维
template<typename Type> using SeVector3 = SeVector<Type, VEC3_XYZ>; 
template<typename Type> using SeVector4 = SeVector<Type, VEC4_XYZW>;

/*************************************************************************
*****************************    SDVector    *****************************
*************************************************************************/

// 向量结构体定义
template<typename Type, int N> struct SeVector
{
	static constexpr int Elements = N; //静态且返回值为常量

	Type m_data[N]; //这里的type可以是多种格式：int\char\...，是m_data中对应元素的格式
	
	 SeVector() 
	{
		for (int i = 0; i < N; ++i)		m_data[i] = 0;
	}

	 SeVector(Type scalar) //如果输入标量，则数组内元素都用该张量赋值
	{
		for (int i = 0; i < N; ++i)		m_data[i] = scalar;
	}

	 //运算符重载：[]访问元素；+-对应元素加减；
	 const Type & operator[](unsigned int i) const { SE_ASSERT(i < Elements); return m_data[i]; }

	 Type & operator[](unsigned int i) { SE_ASSERT(i < Elements); return m_data[i]; }

	 SeVector operator+(const SeVector & a) const 
	{ 
		SeVector ans;
		for (int i = 0; i < N; ++i)		
			ans[i] = m_data[i] + a[i];
		return ans;
	}
	 SeVector operator-(const SeVector & a) const
	{
		SeVector ans;
		for (int i = 0; i < N; ++i)
			ans[i] = m_data[i] - a[i];
		return ans;
	}

	 //数乘
	 SeVector operator*(const Type scalar) const
	{
		SeVector ans;
		for (int i = 0; i < N; ++i)
			ans[i] = m_data[i] * scalar;
		return ans;
	}
	 //数除
	 SeVector operator/(const Type scalar) const //末尾const表明该函数不能修改类内成员
	{
		Type invScalar = Type(1) / scalar;
		SeVector ans;
		for (int i = 0; i < N; ++i)
			ans[i] = m_data[i] * invScalar;
		return ans;
	}

	 void operator+=(const SeVector & a) 
	{
		for (int i = 0; i < N; ++i)
			m_data[i] += a[i];
	}

	 void operator-=(const SeVector & a)
	{
		for (int i = 0; i < N; ++i)
			m_data[i] -= a[i];
	}

	//外部输入一个sevector的数乘运算——友函数！！
	friend  SeVector operator*(Type scalar, const SeVector & a) { return a * scalar; }

	//点乘法
	 Type Dot(const SeVector<Type, N> & rhs) const
	{
		Type sum = 0;
		for (int i = 0; i < N; ++i)
			sum += m_data[i] * rhs.values[i];
		return sum;
	}

	 //求和：每个元素相加
	 Type Sum() const
	{
		Type sum = 0;
		for (int i = 0; i < N; ++i)
			sum += m_data[i];
		return sum;
	}

	 //平方和
	 Type SqrLength()const
	{
		Type sum = 0;
		for (int i = 0; i < N; ++i)
			sum += m_data[i] * m_data[i];
		return sum;
	}

	 //向量长度
	 Type Length()const
	{
		return std::sqrt(SqrLength());
	}
	 //归一化，不修改原始向量
	 SeVector<Type, N> Normalize()const
	{
		return (*this) / Length();
	}
	 //修改原始向量为单位向量
	 SeVector<Type, N> & NormalizeLocal()
	{
		(*this) /= Length();		return *this;
	}
};


//以下对2、3、4维向量结构体进行具体定义：
/*************************************************************************
****************************    SDVector2    *****************************
*************************************************************************/

// 二维向量结构体
template<typename Type> struct SE_VEC2_ALIGN(Type) SeVector<Type, VEC2_XY> //模板特例化——二维情形,前面struct使用了宏定义用于对齐？
{
	static constexpr int Elements = 2;

	using value_type = Type;

	union
	{
		//定义联合体：包含x,y和一个二维数组(可指定类型)；
		struct { Type x, y; };
		struct { Type values[2]; };
	};

	 SeVector() {}
	 SeVector(Type s1, Type s2) : x(s1), y(s2) {}  //根据s1,s2给x,y赋值
	 SeVector(Type scalar) : x(scalar), y(scalar) {} //同样是赋值操作：标量赋值

	 //运算符重载——分为标量输入和向量输入;以及改变原始值和不改变原始值
	 void operator+=(Type scalar) { x += scalar; y += scalar; }
	 void operator-=(Type scalar) { x -= scalar; y -= scalar; }
	 void operator*=(Type scalar) { x *= scalar; y *= scalar; }
	 void operator/=(Type scalar) { x /= scalar; y /= scalar; }

	 void operator+=(const SeVector & a) { x += a.x; y += a.y; }
	 void operator-=(const SeVector & a) { x -= a.x; y -= a.y; }
	 void operator*=(const SeVector & a) { x *= a.x; y *= a.y; }
	 void operator/=(const SeVector & a) { x /= a.x; y /= a.y; }

	 //不改变原始向量的操作
	 SeVector operator+(Type scalar) const { return SeVector(x + scalar, y + scalar); }
	 SeVector operator-(Type scalar) const { return SeVector(x - scalar, y - scalar); }
	 SeVector operator*(Type scalar) const { return SeVector(x * scalar, y * scalar); }
	 SeVector operator/(Type scalar) const { return SeVector(x / scalar, y / scalar); }

	 SeVector operator+(const SeVector & a) const { return SeVector(x + a.x, y + a.y); }
	 SeVector operator-(const SeVector & a) const { return SeVector(x - a.x, y - a.y); }
	 SeVector operator*(const SeVector & a) const { return SeVector(x * a.x, y * a.y); }
	 SeVector operator/(const SeVector & a) const { return SeVector(x / a.x, y / a.y); }

	 //友函数：允许在外部调用，返回其他向量的运算
	friend  SeVector operator+(Type scalar, const SeVector & a) { return SeVector(scalar + a.x, scalar + a.y); }
	friend  SeVector operator-(Type scalar, const SeVector & a) { return SeVector(scalar - a.x, scalar - a.y); }
	friend  SeVector operator*(Type scalar, const SeVector & a) { return SeVector(scalar * a.x, scalar * a.y); }
	friend  SeVector operator/(Type scalar, const SeVector & a) { return SeVector(scalar / a.x, scalar / a.y); }

	//模板函数定义：将vec2/3/4othertype数据格式强制转换为目标格式type后赋值给(x,y)，所以允许输入不同类型的变量进行赋值(进行强制类型转换)
	//结构体内不能定义函数？那这里这个定义操作我不是很明白了，但执行的操作就是支持跨类别的赋值
	template<typename OtherType>  explicit SeVector(const SeVector2<OtherType> & vec2) : x(static_cast<Type>(vec2.x)), y(static_cast<Type>(vec2.y)) {}
	template<typename OtherType>  explicit SeVector(const SeVector3<OtherType> & vec3) : x(static_cast<Type>(vec3.x)), y(static_cast<Type>(vec3.y)) {}
	template<typename OtherType>  explicit SeVector(const SeVector4<OtherType> & vec4) : x(static_cast<Type>(vec4.x)), y(static_cast<Type>(vec4.y)) {}

	//取值运算符[]: 可修改和不可修改
	 const Type & operator[](unsigned int i) const { SE_ASSERT(i < Elements); return values[i]; }
	 Type & operator[](unsigned int i) { SE_ASSERT(i < Elements); return values[i]; }
};

/*************************************************************************
****************************    SDVector3    *****************************
*************************************************************************/
//三元向量类似二元
template<typename Type> struct SeVector<Type, VEC3_XYZ>
{
	static constexpr int Elements = 3;

	using value_type = Type;

	union
	{
		struct { Type x, y, z; };
		struct { Type values[3]; };
	};

	 SeVector() {}
	 SeVector(Type scalar) : x(scalar), y(scalar), z(scalar) {}
	 SeVector(Type s1, Type s2, Type s3) : x(s1), y(s2), z(s3) {}
	 explicit SeVector(Type s1, const SeVector2<Type> & vec2) : x(s1), y(vec2.x), z(vec2.y) {}
	 explicit SeVector(const SeVector2<Type> & vec2, Type s3) : x(vec2.x), y(vec2.y), z(s3) {}

	 void operator+=(const SeVector & a) { x += a.x; y += a.y; z += a.z; }
	 void operator-=(const SeVector & a) { x -= a.x; y -= a.y; z -= a.z; }
	 void operator*=(const SeVector & a) { x *= a.x; y *= a.y; z *= a.z; }
	 void operator/=(const SeVector & a) { x /= a.x; y /= a.y; z /= a.z; }

	 void operator+=(Type scalar) { x += scalar; y += scalar; z += scalar; }
	 void operator-=(Type scalar) { x -= scalar; y -= scalar; z -= scalar; }
	 void operator*=(Type scalar) { x *= scalar; y *= scalar; z *= scalar; }
	 void operator/=(Type scalar) { x /= scalar; y /= scalar; z /= scalar; }

	 SeVector operator+(const SeVector & a) const { return SeVector(x + a.x, y + a.y, z + a.z); }
	 SeVector operator-(const SeVector & a) const { return SeVector(x - a.x, y - a.y, z - a.z); }
	 SeVector operator*(const SeVector & a) const { return SeVector(x * a.x, y * a.y, z * a.z); }
	 SeVector operator/(const SeVector & a) const { return SeVector(x / a.x, y / a.y, z / a.z); }

	 SeVector operator+(Type scalar) const { return SeVector(x + scalar, y + scalar, z + scalar); }
	 SeVector operator-(Type scalar) const { return SeVector(x - scalar, y - scalar, z - scalar); }
	 SeVector operator*(Type scalar) const { return SeVector(x * scalar, y * scalar, z * scalar); }
	 SeVector operator/(Type scalar) const { return SeVector(x / scalar, y / scalar, z / scalar); }

	friend  SeVector operator+(Type scalar, const SeVector & a) { return SeVector(scalar + a.x, scalar + a.y, scalar + a.z); }
	friend  SeVector operator-(Type scalar, const SeVector & a) { return SeVector(scalar - a.x, scalar - a.y, scalar - a.z); }
	friend  SeVector operator*(Type scalar, const SeVector & a) { return SeVector(scalar * a.x, scalar * a.y, scalar * a.z); }
	friend  SeVector operator/(Type scalar, const SeVector & a) { return SeVector(scalar / a.x, scalar / a.y, scalar / a.z); }

	template<typename OtherType>  explicit SeVector(const SeVector2<OtherType> & vec2) : x(static_cast<Type>(vec2.x)), y(static_cast<Type>(vec2.y)), z(static_cast<Type>(0)) {}
	template<typename OtherType>  explicit SeVector(const SeVector3<OtherType> & vec3) : x(static_cast<Type>(vec3.x)), y(static_cast<Type>(vec3.y)), z(static_cast<Type>(vec3.z)) {}
	template<typename OtherType>  explicit SeVector(const SeVector4<OtherType> & vec4) : x(static_cast<Type>(vec4.x)), y(static_cast<Type>(vec4.y)), z(static_cast<Type>(vec4.z)) {}

	 const Type & operator[](unsigned int i) const { SE_ASSERT(i < Elements); return values[i]; }
	 Type & operator[](unsigned int i) { SE_ASSERT(i < Elements); return values[i]; }
};

/*************************************************************************
****************************    SDVector4    *****************************
*************************************************************************/

template<typename Type> struct SE_VEC4_ALIGN(Type) SeVector<Type, VEC4_XYZW>
{
	static constexpr int Elements = 4;

	using value_type = Type;

	union
	{
		struct { Type values[4]; };
		struct { Type x, y, z, w; };
		struct { SeVector3<Type> xyz; }; //多了一项xyz
	};

	 SeVector() {}
	 SeVector(Type scalar) : x(scalar), y(scalar), z(scalar), w(scalar) {}
	 SeVector(Type s1, Type s2, Type s3, Type s4) : x(s1), y(s2), z(s3), w(s4) {}

	 //可以支持通过二元/三元向量对四维向量的赋值操作
	 explicit SeVector(Type s1, Type s2, const SeVector2<Type> & vec2) : x(s1), y(s2), z(vec2.x), w(vec2.y) {}
	 explicit SeVector(Type s1, const SeVector2<Type> & vec2, Type s4) : x(s1), y(vec2.x), z(vec2.y), w(s4) {}
	 explicit SeVector(const SeVector2<Type> & vec2, Type s3, Type s4) : x(vec2.x), y(vec2.y), z(s3), w(s4) {}
	 explicit SeVector(const SeVector2<Type> & vec2a, const SeVector2<Type> & vec2b) : x(vec2a.x), y(vec2a.y), z(vec2b.x), w(vec2b.y) {}
	 explicit SeVector(const SeVector3<Type> & vec3, Type s4) : x(vec3.x), y(vec3.y), z(vec3.z), w(s4) {}
	 explicit SeVector(Type s1, const SeVector3<Type> & vec3) : x(s1), y(vec3.x), z(vec3.y), w(vec3.z) {}
	
	 //可修改结构体内变量的运算：注意什么时候要对最后一维w进行修改；
	 void operator+=(const SeVector3<Type> & a) { x += a.x; y += a.y; z += a.z; }
	 void operator-=(const SeVector3<Type> & a) { x -= a.x; y -= a.y; z -= a.z; }
	 void operator*=(const SeVector3<Type> & a) { x *= a.x; y *= a.y; z *= a.z; }
	 void operator/=(const SeVector3<Type> & a) { x /= a.x; y /= a.y; z /= a.z; }

	 void operator+=(const SeVector & a) { x += a.x; y += a.y; z += a.z; w += a.w; }
	 void operator-=(const SeVector & a) { x -= a.x; y -= a.y; z -= a.z; w -= a.w; }
	 void operator*=(const SeVector & a) { x *= a.x; y *= a.y; z *= a.z; w *= a.w; }
	 void operator/=(const SeVector & a) { x /= a.x; y /= a.y; z /= a.z; w /= a.w; }

	 void operator+=(Type scalar) { x += scalar; y += scalar; z += scalar; w += scalar; }
	 void operator-=(Type scalar) { x -= scalar; y -= scalar; z -= scalar; w -= scalar; }
	 void operator*=(Type scalar) { x *= scalar; y *= scalar; z *= scalar; w *= scalar; }
	 void operator/=(Type scalar) { x /= scalar; y /= scalar; z /= scalar; w /= scalar; }

	 //不可修改结构体内元素的操作
	 SeVector operator+(const SeVector & a) const { return SeVector(x + a.x, y + a.y, z + a.z, w + a.w); }
	 SeVector operator-(const SeVector & a) const { return SeVector(x - a.x, y - a.y, z - a.z, w - a.w); }
	 SeVector operator*(const SeVector & a) const { return SeVector(x * a.x, y * a.y, z * a.z, w * a.w); }
	 SeVector operator/(const SeVector & a) const { return SeVector(x / a.x, y / a.y, z / a.z, w / a.w); }

	 SeVector operator+(const SeVector3<Type> & a) const { return SeVector(x + a.x, y + a.y, z + a.z, w); }
	 SeVector operator-(const SeVector3<Type> & a) const { return SeVector(x - a.x, y - a.y, z - a.z, w); }
	 SeVector operator*(const SeVector3<Type> & a) const { return SeVector(x * a.x, y * a.y, z * a.z, w); }
	 SeVector operator/(const SeVector3<Type> & a) const { return SeVector(x / a.x, y / a.y, z / a.z, w); }

	 SeVector operator+(Type scalar) const { return SeVector(x + scalar, y + scalar, z + scalar, w + scalar); }
	 SeVector operator-(Type scalar) const { return SeVector(x - scalar, y - scalar, z - scalar, w - scalar); }
	 SeVector operator*(Type scalar) const { return SeVector(x * scalar, y * scalar, z * scalar, w * scalar); }
	 SeVector operator/(Type scalar) const { return SeVector(x / scalar, y / scalar, z / scalar, w / scalar); }

	friend  SeVector operator+(Type scalar, const SeVector & a) { return SeVector(scalar + a.x, scalar + a.y, scalar + a.z, scalar + a.w); }
	friend  SeVector operator-(Type scalar, const SeVector & a) { return SeVector(scalar - a.x, scalar - a.y, scalar - a.z, scalar - a.w); }
	friend  SeVector operator*(Type scalar, const SeVector & a) { return SeVector(scalar * a.x, scalar * a.y, scalar * a.z, scalar * a.w); }
	friend  SeVector operator/(Type scalar, const SeVector & a) { return SeVector(scalar / a.x, scalar / a.y, scalar / a.z, scalar / a.w); }

	template<typename OtherType>  explicit SeVector(const SeVector2<OtherType> & vec2) : x(static_cast<Type>(vec2.x)), y(static_cast<Type>(vec2.y)), z(static_cast<Type>(0)), w(static_cast<Type>(0)) {}
	template<typename OtherType>  explicit SeVector(const SeVector3<OtherType> & vec3) : x(static_cast<Type>(vec3.x)), y(static_cast<Type>(vec3.y)), z(static_cast<Type>(vec3.z)), w(static_cast<Type>(0)) {}
	template<typename OtherType>  explicit SeVector(const SeVector4<OtherType> & vec4) : x(static_cast<Type>(vec4.x)), y(static_cast<Type>(vec4.y)), z(static_cast<Type>(vec4.z)), w(static_cast<Type>(vec4.w)) {}

	 const Type & operator[](unsigned int i) const { SE_ASSERT(i < Elements); return values[i]; }
	 Type & operator[](unsigned int i) { SE_ASSERT(i < Elements); return values[i]; }
};

/*************************************************************************
****************************    Operators    *****************************
*************************************************************************/

//运算符定义：逻辑运算符等
template<typename Type, int N> SE_INLINE bool operator==(const SeVector<Type, N> & a, const SeVector<Type, N> & b)
{
	for (int i = 0; i < SeVector<Type, N>::Elements; ++i)		{ if (a[i] != b[i])	return false; }			return true;
}
template<typename Type, int N> SE_INLINE bool operator!=(const SeVector<Type, N> & a, const SeVector<Type, N> & b)
{
	for (int i = 0; i < SeVector<Type, N>::Elements; ++i)		{ if (a[i] != b[i])	return true; }			return false;
}
template<typename Type, int N> SE_INLINE SeVector<Type, N> operator-(SeVector<Type, N> a)
{
	for (int i = 0; i < SeVector<Type, N>::Elements; ++i)		a[i] = -a[i];			return a;
}
template<typename Type, int N> SE_INLINE SeVector<Type, N> operator!(SeVector<Type, N> a)
{
	for (int i = 0; i < SeVector<Type, N>::Elements; ++i)		a[i] = !a[i];			return a;
}
template<typename Type, int N> SE_INLINE SeVector<Type, N> operator~(SeVector<Type, N> a)
{
	for (int i = 0; i < SeVector<Type, N>::Elements; ++i)		a[i] = ~a[i];			return a;
}


/*************************************************************************
***************************    Type defines    ***************************
*************************************************************************/
//变量类型定义：左边的类型统一用右边的名字来命名
// bool
typedef SeVector2<bool>				Bool2;
// char
typedef SeVector2<char>				Char2;
typedef SeVector3<char>				Char3;
typedef SeVector4<char>				Char4;
typedef SeVector<char, 5>			Char5;
typedef SeVector<char, 6>			Char6;
// unsign char
typedef SeVector2<unsigned char>	UChar2;
typedef SeVector3<unsigned char>	UChar3;
typedef SeVector4<unsigned char>	UChar4;
typedef SeVector<unsigned char, 5>	UChar5;
typedef SeVector<unsigned char, 6>	UChar6;
// short
typedef SeVector2<short>			Short2;
typedef SeVector3<short>			Short3;
typedef SeVector4<short>			Short4;
typedef SeVector<short, 5>			Short5;
typedef SeVector<short, 6>			Short6;
// unsigned short
typedef SeVector2<unsigned short>	UShort2;
typedef SeVector3<unsigned short>	UShort3;
typedef SeVector4<unsigned short>	UShort4;
typedef SeVector<unsigned short, 5> UShort5;
typedef SeVector<unsigned short, 6> UShort6;
// int
typedef SeVector2<int>				Int2;
typedef SeVector3<int>				Int3;
typedef SeVector4<int>				Int4;
typedef SeVector<int, 5>			Int5;
typedef SeVector<int, 6>			Int6;
typedef SeVector<int, 8>			Int8;
typedef SeVector<int, 12>			Int12;
typedef SeVector<int, 30>			Int30;
// unsigned int
typedef SeVector2<unsigned int>		UInt2;
typedef SeVector3<unsigned int>		UInt3;
typedef SeVector4<unsigned int>		UInt4;
typedef SeVector<unsigned int, 5>	UInt5;
typedef SeVector<unsigned int, 6>	UInt6;
// float
typedef SeVector2<float>			Float2;
typedef SeVector3<float>			Float3;
typedef SeVector4<float>			Float4;
typedef SeVector<float, 5>			Float5;
typedef SeVector<float, 6>			Float6;
typedef SeVector<float, 9>			Float9;
typedef SeVector<float, 12>			Float12;
// double
typedef SeVector2<double>			Double2;
typedef SeVector3<double>			Double3;
typedef SeVector4<double>			Double4;
typedef SeVector<double, 5>			Double5;
typedef SeVector<double, 6>			Double6;
// size_t
typedef SeVector2<size_t>			Size2;
typedef SeVector3<size_t>			Size3;
typedef SeVector4<size_t>			Size4;
typedef SeVector<size_t, 5>			Size5;
typedef SeVector<size_t, 6>			Size6;
/*
// SIMD-float3
typedef SeVector<float, SIMD_VEC3>	Simd3f;*/

typedef SeVector3<short>			Half3;	//just for compile reason,don't use it 

SE_NAMESPACE_END