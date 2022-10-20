#pragma once
#include "edge_face.h"

namespace my
{
	void gen_edge(const int numvert, const int row_num, const int col_num, SE::Int4* m_edges, SE::SeVec3fSimd* m_positions)
	{
		//SE::Int4* m_edges = new SE::Int4[numedge];
		int ecount = 0;
		for (int i = 0; i < numvert; i++)
		{
			if (m_positions[i].x < (float)(row_num-1)*0.25 && m_positions[i].y < (float)(col_num - 1) * 0.25 )
			{
				SE::Int4 temp1 = { i, i + 1, 0, 0 };
				SE::Int4 temp2 = { i, i + row_num, 0, 0 };
				SE::Int4 temp3 = { i, i + row_num+1, 0, 0 };
				m_edges[ecount] = temp1;
				m_edges[ecount + 1] = temp2;
				m_edges[ecount + 2] = temp3;
				ecount += 3;
			}
			else
			{
				if (m_positions[i].x == (float)(row_num - 1) * 0.25 && m_positions[i].y < (float)(col_num - 1) * 0.25 )
				{
					SE::Int4 temp1 = { i, i + 1, 0, 0 };
					m_edges[ecount] = temp1;
					ecount += 1;
				}
				if (m_positions[i].y == (float)(col_num - 1) * 0.25 && m_positions[i].x < (float)(row_num - 1) * 0.25 )
				{
					SE::Int4 temp1 = { i, i + row_num, 0, 0 };
					m_edges[ecount] = temp1;
					ecount += 1;
				}
			}
		}
	}

	void gen_face(const int numvert, const int row_num, const int col_num, SE::Int4* m_faces, SE::SeVec3fSimd* m_positions)
	{
		//SE::Int4* m_edges = new SE::Int4[numedge];
		int fcount = 0;
		for (int i = 0; i < numvert; i++)
		{
			if (m_positions[i].x == (float)(row_num - 1) * 0.25 || m_positions[i].y == (float)(col_num - 1) * 0.25) { continue; }
			SE::Int4 temp2 = { i, i + row_num, i + col_num+1 , 0 };
			SE::Int4 temp1 = { i, i + 1, i + col_num + 1, 0 };
			m_faces[fcount] = temp1;
			m_faces[fcount + 1] = temp2;
			fcount += 2;
		}
	}

	int neibor_start(const int numvert, const int row_num, const int col_num, int* start)
	{
		start[0]=0;
		for (unsigned int i = 1; i < numvert + 1; i++)
		{
			int temp = start[i - 1];
			//first row
			if ((i - 1) / (row_num) == 0) 
			{
				if (i == 1)
				{
					temp += 3;
					start[i]=temp;
					continue;
				}
				if (i == row_num)
				{
					temp += 2;
					start[i]=temp;
					continue;
				}
				temp += 4;
				start[i]=temp;
				continue;
			}

			if ((i - 1) / (row_num) < (col_num-1)) //2->(last_row-1)
			{
				if ((i - 1) % (row_num) == 0)
				{
					temp += 4;
					start[i]=temp;
					continue;
				}
				if ((i - 1) % ((row_num)) == ((row_num)-1))
				{
					temp += 4;
					start[i]=temp;
					continue;
				}
				temp += 6;
				start[i]=temp;
				continue;
			}
			if ((i - 1) / (row_num) == (col_num-1))
			{
				if ((i - 1) % (row_num) == 0)
				{
					temp += 2;
					start[i]=temp;
					continue;
				}
				if ((i - 1) % (row_num) == (row_num-1))
				{
					temp += 3;
					start[i]=temp;
					continue;
				}
				temp += 4;
				start[i]=temp;
			}
		}
		return start[numvert];
	}
	void neibor_idx(const int numvert, const int row_num, const int col_num, int *idx)
	{
		int temp = 0;
		for (int i = 1; i < numvert + 1; i++)
		{
			//first row
			if ((i - 1) / (row_num) == 0)
			{
				if (i == 1)
				{
					idx[temp] = i;
					idx[temp+1] = i + row_num-1;
					idx[temp+2] = i + row_num;
					temp += 3;
					continue;
				}
				if (i == row_num)
				{
					idx[temp] = i - 2;
					idx[temp+1] = i + row_num-1;
					temp += 2;
					continue;
				}
				idx[temp] = i - 2;
				idx[temp + 1] = i ;
				idx[temp + 2] = i+ row_num - 1;
				idx[temp + 3] = i+ row_num;
				temp += 4;
				continue;
			}
			if ((i - 1) / row_num < (col_num - 1))
			{
				if ((i - 1) % row_num == 0)
				{
					idx[temp] = i - row_num - 1;
					idx[temp + 1] = i;
					idx[temp + 2] = i + row_num - 1;
					idx[temp + 3] = i + row_num;
					temp += 4;
					continue;
				}
				if ((i - 1) % row_num == row_num - 1)
				{
					idx[temp] = i - row_num - 2;
					idx[temp + 1] = i - row_num - 1;
					idx[temp + 2] = i - 2;
					idx[temp + 3] = i + row_num - 1;
					temp += 4;
					continue;
				}
				idx[temp] = i - row_num - 2;
				idx[temp + 1] = i - row_num - 1;
				idx[temp + 2] = i - 2;
				idx[temp + 3] = i;
				idx[temp + 4] = i + row_num - 1;
				idx[temp + 5] = i + row_num;
				temp += 6;
				continue;
			}
			if ((i - 1) / row_num == col_num - 1)
			{
				if ((i - 1) % row_num == 0)
				{
					idx[temp] = i - row_num - 1;
					idx[temp + 1] = i;
					temp += 2;
					continue;
				}
				if ((i - 1) % row_num == row_num - 1)
				{
					idx[temp] = i - row_num - 2;
					idx[temp + 1] = i - row_num - 1;
					idx[temp + 2] = i - 2;
					temp += 3;	
					continue;
				}
				idx[temp] = i - row_num - 2;
				idx[temp + 1] = i - row_num - 1;
				idx[temp + 2] = i - 2;
				idx[temp + 3] = i;
				temp += 4;
			}
		}	
	}
}