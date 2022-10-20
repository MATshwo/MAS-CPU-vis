#pragma once
#include "SeVector.h"
#include "SeVectorSimd.h"
#include "SeMatrix.h"
#include "SeCsr.h"
#include "SeCollisionElements.h"
namespace my 
{
	void gen_edge(const int numvert, const int row_num, const int col_num, SE::Int4* m_edges, SE::SeVec3fSimd* m_positions);
	void gen_face(const int numvert, const int row_num, const int col_num, SE::Int4* m_faces, SE::SeVec3fSimd* m_positions);
	int neibor_start(const int numvert, const int row_num, const int col_num, int* start);
	void neibor_idx(const int numvert, const int row_num, const int col_num, int* idx);
}