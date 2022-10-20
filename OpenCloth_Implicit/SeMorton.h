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

SE_NAMESPACE_BEGIN

/*************************************************************************
****************************    SeMorton64    ****************************
*************************************************************************/

/**
 *	@brief	64-bit Morton-Code object. 64位的莫顿编码
 */
// LL是long long int , 1ll就是长整型的数值1
//代码里处处使用value_type(),这点我还存在疑惑？吃完饭回来再了解下！
class SE_ALIGN(8) SeMorton64
{

public:
	//64位无符号整型：最大可表示2^64-1，按理论来说，莫顿码是63位的无符号整型/64的整型数据(多一位表示符号)
	using value_type = unsigned long long; //这里using关键字的用作重命名

	//!	@brief	Default constructor.
	 SeMorton64() {} //构造函数-无输入，与类名相同

	//!	@brief	Convert to value_type.
	//这里operator不是用于运算符重载，而是隐式类型转换，将类输出由其他格式最终转换为value_type类型并返回；const关键字表示常量成员函数
	 operator value_type() const { return m_Value; }  //隐式类型转换？

	//!	@brief	Constructed by a given value. 根据给定值构造
	 SeMorton64(value_type _Value) : m_Value(_Value) {} //构造函数重载-有参数输入的版本

	//!	@brief	Encoded by the given 3d point located within the unit cube (controlable precision).由三维坐标求其莫顿编码的模板函数
	template<unsigned int precision>  void Encode(float x, float y, float z, value_type lastBits) 
	{
		// 可指定precision位(坐标二进制位：要求<21 )的莫顿编码函数
		static_assert(precision <= 21, "The highest precision for 64-bit Morton code is 21."); 

		x = Math::Clamp(x * (1ll << precision), 0.0f, (1ll << precision) - 1.0f); //1ll << precision，将长整型数据1左移precision位，clamp先比大后比小，最大不超过2^21-1,长度就是21.
		y = Math::Clamp(y * (1ll << precision), 0.0f, (1ll << precision) - 1.0f);
		z = Math::Clamp(z * (1ll << precision), 0.0f, (1ll << precision) - 1.0f);

		value_type xx = SeMorton64::ExpandBits(static_cast<value_type>(x)) << (66 - 3 * precision); //不足64位补齐的操作
		value_type yy = SeMorton64::ExpandBits(static_cast<value_type>(y)) << (65 - 3 * precision);
		value_type zz = SeMorton64::ExpandBits(static_cast<value_type>(z)) << (64 - 3 * precision);

		constexpr value_type bitMask = ~value_type(0) >> (3 * precision);

		m_Value = xx + yy + zz + (lastBits & bitMask);
	}

	//!	@brief	Encoded by the given 3d point located within the unit cube (full precision). 默认21位版本
	 void Encode(float x, float y, float z) 
	{
		//这里应该是预处理后的坐标，默认是[0,1]之间的小数，所以这里的操作是将小数后的位进行平移到小数点之前
		//即数据缩放：[0,1]-> [0,2^21-1];括号里三个值，第一个表示平移(即缩放操作),后两个是为了限制数据的范围，确保数据在这个范围内。
		x = Math::Clamp(x * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);  //这里还是浮点型
		y = Math::Clamp(y * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);
		z = Math::Clamp(z * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);

		//这里就已经强制转换为ll类型了：static_cast<type>(x) 将x强制转类型
		value_type xx = SeMorton64::ExpandBits(static_cast<value_type>(x)); //将原始二进制扩充：每两位和首位前增加两个0
		value_type yy = SeMorton64::ExpandBits(static_cast<value_type>(y));
		value_type zz = SeMorton64::ExpandBits(static_cast<value_type>(z));

		m_Value = (xx << 2) + (yy << 1) + zz; //返回结果再左移2位，y左移1位，z不移动（奇怪这个和网上不太一样，网上z在最左边）
	}

private:

	/**
	 *	@brief	Expand bits by inserting two zeros after each bit. //扩展位：在每两位之间插入两位0位,21位变成了63位，注意首项之前还需要再补两个0.
	 *	@e.g.	0000 0000 1111  ->  0010 0100 1001  4位变12位
	 */
	static  value_type ExpandBits(value_type bits)
	{
		//进行五组运算：每一组都是由平移+按位与+按位或的运算构成！//按道理进行四次就可以欸，不知道为啥是五次，好像第一步可以去掉，判断是否32，默认输入32位
		bits = (bits | (bits << 32)) & 0xFFFF00000000FFFFu; 
		bits = (bits | (bits << 16)) & 0x00FF0000FF0000FFu;
		bits = (bits | (bits <<  8)) & 0xF00F00F00F00F00Fu;
		bits = (bits | (bits <<  4)) & 0x30C30C30C30C30C3u;
		return (bits | (bits <<  2)) & 0x9249249249249249u;
	}
private:

	value_type	m_Value;
};

SE_NAMESPACE_END