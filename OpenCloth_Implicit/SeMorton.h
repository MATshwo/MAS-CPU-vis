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
 *	@brief	64-bit Morton-Code object. 64λ��Ī�ٱ���
 */
// LL��long long int , 1ll���ǳ����͵���ֵ1
//�����ﴦ��ʹ��value_type(),����һ������ɻ󣿳��극�������˽��£�
class SE_ALIGN(8) SeMorton64
{

public:
	//64λ�޷������ͣ����ɱ�ʾ2^64-1����������˵��Ī������63λ���޷�������/64����������(��һλ��ʾ����)
	using value_type = unsigned long long; //����using�ؼ��ֵ�����������

	//!	@brief	Default constructor.
	 SeMorton64() {} //���캯��-�����룬��������ͬ

	//!	@brief	Convert to value_type.
	//����operator����������������أ�������ʽ����ת�������������������ʽ����ת��Ϊvalue_type���Ͳ����أ�const�ؼ��ֱ�ʾ������Ա����
	 operator value_type() const { return m_Value; }  //��ʽ����ת����

	//!	@brief	Constructed by a given value. ���ݸ���ֵ����
	 SeMorton64(value_type _Value) : m_Value(_Value) {} //���캯������-�в�������İ汾

	//!	@brief	Encoded by the given 3d point located within the unit cube (controlable precision).����ά��������Ī�ٱ����ģ�庯��
	template<unsigned int precision>  void Encode(float x, float y, float z, value_type lastBits) 
	{
		// ��ָ��precisionλ(���������λ��Ҫ��<21 )��Ī�ٱ��뺯��
		static_assert(precision <= 21, "The highest precision for 64-bit Morton code is 21."); 

		x = Math::Clamp(x * (1ll << precision), 0.0f, (1ll << precision) - 1.0f); //1ll << precision��������������1����precisionλ��clamp�ȱȴ���С����󲻳���2^21-1,���Ⱦ���21.
		y = Math::Clamp(y * (1ll << precision), 0.0f, (1ll << precision) - 1.0f);
		z = Math::Clamp(z * (1ll << precision), 0.0f, (1ll << precision) - 1.0f);

		value_type xx = SeMorton64::ExpandBits(static_cast<value_type>(x)) << (66 - 3 * precision); //����64λ����Ĳ���
		value_type yy = SeMorton64::ExpandBits(static_cast<value_type>(y)) << (65 - 3 * precision);
		value_type zz = SeMorton64::ExpandBits(static_cast<value_type>(z)) << (64 - 3 * precision);

		constexpr value_type bitMask = ~value_type(0) >> (3 * precision);

		m_Value = xx + yy + zz + (lastBits & bitMask);
	}

	//!	@brief	Encoded by the given 3d point located within the unit cube (full precision). Ĭ��21λ�汾
	 void Encode(float x, float y, float z) 
	{
		//����Ӧ����Ԥ���������꣬Ĭ����[0,1]֮���С������������Ĳ����ǽ�С�����λ����ƽ�Ƶ�С����֮ǰ
		//���������ţ�[0,1]-> [0,2^21-1];����������ֵ����һ����ʾƽ��(�����Ų���),��������Ϊ���������ݵķ�Χ��ȷ�������������Χ�ڡ�
		x = Math::Clamp(x * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);  //���ﻹ�Ǹ�����
		y = Math::Clamp(y * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);
		z = Math::Clamp(z * (1ll << 21), 0.0f, (1ll << 21) - 1.0f);

		//������Ѿ�ǿ��ת��Ϊll�����ˣ�static_cast<type>(x) ��xǿ��ת����
		value_type xx = SeMorton64::ExpandBits(static_cast<value_type>(x)); //��ԭʼ���������䣺ÿ��λ����λǰ��������0
		value_type yy = SeMorton64::ExpandBits(static_cast<value_type>(y));
		value_type zz = SeMorton64::ExpandBits(static_cast<value_type>(z));

		m_Value = (xx << 2) + (yy << 1) + zz; //���ؽ��������2λ��y����1λ��z���ƶ��������������ϲ�̫һ��������z������ߣ�
	}

private:

	/**
	 *	@brief	Expand bits by inserting two zeros after each bit. //��չλ����ÿ��λ֮�������λ0λ,21λ�����63λ��ע������֮ǰ����Ҫ�ٲ�����0.
	 *	@e.g.	0000 0000 1111  ->  0010 0100 1001  4λ��12λ
	 */
	static  value_type ExpandBits(value_type bits)
	{
		//�����������㣺ÿһ�鶼����ƽ��+��λ��+��λ������㹹�ɣ�//����������ĴξͿ��ԚG����֪��Ϊɶ����Σ������һ������ȥ�����ж��Ƿ�32��Ĭ������32λ
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