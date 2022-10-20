//#include<iostream>
//#include"SeSchwarzPreconditioner.h"
//#include<emmintrin.h>

/*int main()
{
	//指定SE命名空间,否则无法检测
	//SE::SeArray2D<int> B(2, 2, 1);
	//std::cout << B.Rows() << std::endl;

	//Loading data
	//1、AllocatePrecoditioner input
	const int numvert = 64; //8*8
	const int numedge = 161;  //7*7*3+（8-1）*2   除最后一行&最后一列，其余顶点各自导出三条边，最后一行/列导出行/列顶点数-1条边
	const int numface = 98; //7*7*2 除最后一行&最后一列，其余顶点各自导出两个三角形



	//顶点(x,y,z,0),长度=m_numvert=4,指定四个顶点坐标——考虑增加到64个顶点
	//const SE::SeVec3fSimd* m_positions = new SE::SeVec3fSimd[65]{ {0.0f,0.0f,0.0f},{0.0f,1.0f,0.0f},{1.0f,0.0f,0.0f},{1.0f,1.0f,0.0f} }; //添加const之后就不能随意变动了，最好不要加
	//初始化顶点坐标
	SE::SeVec3fSimd* m_positions = new SE::SeVec3fSimd[numvert];
	//还是赋值的问题
	int vcount = 0;
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{

			SE::SeVec3fSimd temp((float)i / 4, (float)j / 4, 0, 0);
			m_positions[vcount] = temp;
			vcount++;
		}
	}

	//边坐标初始化(e1,e2,0,0),长度=m_numedge
	//const SE::Int4* m_edges = new SE::Int4[numedge]{ {0,1,0,0},{0,2,0,0},{0,3,0,0},{1,3,0,0},{2,3,0,0} };
	SE::Int4* m_edges = new SE::Int4[numedge];
	int ecount = 0;
	for (int i = 0; i < numvert; i++)
	{
		if (m_positions[i].x < 1.75 && m_positions[i].y < 1.75)
		{
			SE::Int4 temp1 = { i, i + 1, 0, 0 };
			SE::Int4 temp2 = { i, i + 8, 0, 0 };
			SE::Int4 temp3 = { i, i + 9, 0, 0 };
			m_edges[ecount] = temp1;
			m_edges[ecount + 1] = temp2;
			m_edges[ecount + 2] = temp3;
			ecount += 3;
		}
		else
		{
			if (m_positions[i].x == 1.75 && m_positions[i].y < 1.75)
			{
				SE::Int4 temp1 = { i, i + 1, 0, 0 };
				m_edges[ecount] = temp1;
				ecount += 1;
			}
			if (m_positions[i].y == 1.75 && m_positions[i].x < 1.75)
			{
				SE::Int4 temp1 = { i, i + 8, 0, 0 };
				m_edges[ecount] = temp1;
				ecount += 1;
			}
		}
	}

	//std::cout << m_edges[0][0] << std::endl;
	//面坐标初始化(x1,x2,x3,0),长度=m_numface
	//const SE::Int4* m_faces = new SE::Int4[2]{ {0,1,3,0},{0,2,3,0} };
	SE::Int4* m_faces = new SE::Int4[numface];
	int fcount = 0;
	for (int i = 0; i < numvert; i++)
	{
		if (m_positions[i].x == 1.75 || m_positions[i].y == 1.75) { continue; }
		SE::Int4 temp2 = { i, i + 8, i + 9, 0 };
		SE::Int4 temp1 = { i, i + 1, i + 9, 0 };
		m_faces[fcount] = temp1;
		m_faces[fcount + 1] = temp2;
		fcount += 2;
	}

	//neighbor信息
	//const SE::SeCsr<int>* m_neighbours = new SE::SeCsr<int>[4];
	//std::vector<int> start = { 0, 2, 4, 6 }; //run failed 长度不对所以报错 = numvert+1
	std::vector<int> start; //run auccess
	std::vector<int> idx;
	start.push_back(0);

	for (unsigned int i = 1; i < 64 + 1; i++)
	{
		int temp = start[i - 1];
		//first row
		if ((i - 1) / 8 == 0)
		{
			if (i == 1)
			{
				temp += 3;
				start.push_back(temp);
				idx.push_back(i);
				idx.push_back(i + 7);
				idx.push_back(i + 8);
				continue;
			}
			if (i == 8)
			{
				temp += 2;
				start.push_back(temp);
				idx.push_back(i - 2);
				idx.push_back(i + 7);
				continue;
			}
			temp += 4;
			start.push_back(temp);
			idx.push_back(i - 2);
			idx.push_back(i);
			idx.push_back(i + 7);
			idx.push_back(i + 8);

			continue;
		}
		if ((i - 1) / 8 < 7)
		{
			if ((i - 1) % 8 == 0)
			{
				temp += 4;
				start.push_back(temp);
				idx.push_back(i - 9);
				idx.push_back(i);
				idx.push_back(i + 7);
				idx.push_back(i + 8);
				continue;
			}
			if ((i - 1) % 8 == 7)
			{
				temp += 4;
				start.push_back(temp);
				idx.push_back(i - 10);
				idx.push_back(i - 9);
				idx.push_back(i - 2);
				idx.push_back(i + 7);
				continue;
			}
			temp += 6;
			start.push_back(temp);
			idx.push_back(i - 10);
			idx.push_back(i - 9);
			idx.push_back(i - 2);
			idx.push_back(i);
			idx.push_back(i + 7);
			idx.push_back(i + 8);

			continue;
		}
		if ((i - 1) / 8 == 7)
		{
			if ((i - 1) % 8 == 0)
			{
				temp += 2;
				start.push_back(temp);
				idx.push_back(i - 9);
				idx.push_back(i);
				continue;
			}
			if ((i - 1) % 8 == 7)
			{
				temp += 3;
				start.push_back(temp);
				idx.push_back(i - 10);
				idx.push_back(i - 9);
				idx.push_back(i - 2);
				continue;
			}
			temp += 4;
			start.push_back(temp);
			idx.push_back(i - 10);
			idx.push_back(i - 9);
			idx.push_back(i - 2);
			idx.push_back(i);
		}
	}
	std::vector<int> values(start.back(), 1);  //长度=start的最后一个元素：全1向量,存储所有的邻接信息  : start.back() = idx.size()

	const SE::SeCsr<int> m_neighbours_temp(start, idx, values); //const SE::SeCsr<int> m_neighbours(start,idx,values)
	const SE::SeCsr<int>* m_neighbours = &m_neighbours_temp; //指针变量赋值的方法: 创建普通变量后取地址赋值;或者使用new整体赋值
	//std::cout << m_neighbours->Size() << std::endl; //指针变量用->,普通变量用.

	//2、PreparePreconditioner input
	//length = numverts
	//const SE::SeMatrix3f* diagonal = new SE::SeMatrix3f[numvert]{ SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };
	SE::SeMatrix3f* diagonal = new SE::SeMatrix3f[numvert]{ SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };

	//csrrange类似start,存储每个质点首个offdiagnoal在列表csrrange[]中的index  
	int* csrRanges = new int[numvert + 1]; //表示第i个顶点的前i-1个顶点一共有多少邻居，按理说应该=start
	int nn = 0;
	for (auto i = start.begin(); i != start.end(); i++) { csrRanges[nn] = (int)*i; nn++; }

	//std::cout << csrRanges[0] << std::endl;
	//存储质点之间的连接信息,length = sum_{i}(neighbor[i]-1):所有顶点连接的邻居个数的总和
	//SE::SeMatrix3f* csrOffDiagonals = new SE::SeMatrix3f[start.back()]{ SE::SeMatrix3f::Identity(3),SE::SeMatrix3f::Identity(1),SE::SeMatrix3f::Identity(2),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };
	SE::SeMatrix3f* csrOffDiagonals = new SE::SeMatrix3f[start.back()];
	
	//matrix3f赋值运算：两种方式，如果需要逐个指定元素，用float list[9]赋值，如果只是主对角线赋值：SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(a+1);
	//for (int a = 0; a < start.back(); a++)
	//{
	//	float r[9] = { 0,1,2,3,4,5,6,7,8 };
	//	SE::Float3 temp(1,2,3);
	//	SE::SeMatrix3f r1(r);//注意这里是按列赋值，不是按行
	//	SE::SeMatrix3f r2 = SE::SeMatrix3f::Identity(1);
	//	SE::SeMatrix3f r3; 		r3.SetDiagonal(temp);
	//}

	for (int a = 0; a < start.back(); a++)
	{
		SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(0); //非对角元素指定不同导致结果的不同
		csrOffDiagonals[a] = nows;
	}
	//参考论文PSCC: 此处的集合需要通过碰撞算法检测得到,表示下一时刻检测到碰撞的所有情形
	const SE::EfSet* efSets = {}; //先不考虑ef碰撞
	const SE::EeSet* eeSets = {};
	const SE::VfSet* vfSets = {};
	unsigned int* efCounts = new unsigned int[numedge + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; //整数列表初始化为0
	unsigned int* eeCounts = new unsigned int[numedge + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	unsigned int* vfCounts = new unsigned int[numvert + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	//3、Preconditioning input
	//std::vector<SE::SeVec3fSimd> z;
	SE::SeVec3fSimd* z = new SE::SeVec3fSimd[numvert]; //初始为空,用于存储计算后的结果
	
	//residual 右侧常数项
	SE::SeVec3fSimd* residual = new SE::SeVec3fSimd[numvert]; //之后有现成结果赋值再使用const定义，这里先删掉了
	for (int a = 0; a < numvert; a++)
	{
		SE::SeVec3fSimd nows(0.2, 0.1, 0.3, 0);
		residual[a] = nows;
	}
	int dim = 1;
	//4、prepare preconditioner 
	//SE::SeSchwarzPreconditioner* sh = nullptr; 
	SE::SeSchwarzPreconditioner sh;
	sh.m_positions = m_positions;
	sh.m_edges = m_edges;
	sh.m_faces = m_faces;
	sh.m_neighbours = m_neighbours; //若neighbor为普通变量,此处赋值添加&

	sh.AllocatePrecoditioner(numvert, numedge, numface);
	std::cout << "AllocatePrecoditioner run success!" << std::endl;
	//
	///PCG迭代时调用：放在PCG循环里，每次循环得到系统矩阵
	sh.PreparePreconditioner(diagonal, csrOffDiagonals, csrRanges, efSets, eeSets, vfSets, efCounts, eeCounts, vfCounts);
	std::cout << "PreparePreconditioner run success!" << std::endl;
	//
	sh.Preconditioning(z, residual, dim);
	std::cout << "Preconditioning run success!" << std::endl;

	////show result
	std::cout << "Final Z:\n" << std::endl;
	//for (int k = 0; k < numvert; k++)
	//{
	//	printf("%6f %6f %6f\n", z[k].x, z[k].y, z[k].z);
	//}
	return 0;
}
*/