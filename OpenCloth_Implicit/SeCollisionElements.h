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

#include "SeVectorSimd.h"

SE_NAMESPACE_BEGIN

struct EfSet
{
	int		m_eId=0; //edge index
	int		m_fId=0; //face index
	float	stiff=0; //stiff
	Float3	m_bary = {0,0,0};			// (x, 1-x) / (y, z, 1-y-z): barycentric weight of the intersection point 交点的重心权重
	SeVec3fSimd m_normal = {0,0,0};	// repulsion direction 斥力方向
};

struct VfSet
{
	int		m_vId;
	int		m_fId;
	float	stiff;
	Float2	m_bary;			// (x, y, 1-x-y): barycentric weight of the vertex 顶点重心权重
	SeVec3fSimd m_normal;	// repulsion direction 
};

struct EeSet
{
	int		m_eId0;
	int		m_eId1;
	float	stiff;
	Float2	m_bary;			// (x, 1-x) / (y, 1-y): barycentric weight of the two closest points 两个最近点的重心权重
	SeVec3fSimd m_normal;	// repulsion direction
};

struct Stencil  //模板
{
	int verextNumPerStencil; //每个stencil的顶点数目
	int vertexNumOfFirstPrimitive;  //初始的顶点个数

	int		index[5]; //长度为5的列表：原因,最大长度=5
	float	weight[5]; //权重：每个index贡献的权重
	float	stiff; //刚度系数
	SeVec3fSimd direction; //法线方向
};

SE_NAMESPACE_END