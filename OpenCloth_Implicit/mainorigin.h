//#include<iostream>
//#include"SeSchwarzPreconditioner.h"
//#include<emmintrin.h>

/*int main()
{
	//ָ��SE�����ռ�,�����޷����
	//SE::SeArray2D<int> B(2, 2, 1);
	//std::cout << B.Rows() << std::endl;

	//Loading data
	//1��AllocatePrecoditioner input
	const int numvert = 64; //8*8
	const int numedge = 161;  //7*7*3+��8-1��*2   �����һ��&���һ�У����ඥ����Ե��������ߣ����һ��/�е�����/�ж�����-1����
	const int numface = 98; //7*7*2 �����һ��&���һ�У����ඥ����Ե�������������



	//����(x,y,z,0),����=m_numvert=4,ָ���ĸ��������ꡪ���������ӵ�64������
	//const SE::SeVec3fSimd* m_positions = new SE::SeVec3fSimd[65]{ {0.0f,0.0f,0.0f},{0.0f,1.0f,0.0f},{1.0f,0.0f,0.0f},{1.0f,1.0f,0.0f} }; //���const֮��Ͳ�������䶯�ˣ���ò�Ҫ��
	//��ʼ����������
	SE::SeVec3fSimd* m_positions = new SE::SeVec3fSimd[numvert];
	//���Ǹ�ֵ������
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

	//�������ʼ��(e1,e2,0,0),����=m_numedge
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
	//�������ʼ��(x1,x2,x3,0),����=m_numface
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

	//neighbor��Ϣ
	//const SE::SeCsr<int>* m_neighbours = new SE::SeCsr<int>[4];
	//std::vector<int> start = { 0, 2, 4, 6 }; //run failed ���Ȳ������Ա��� = numvert+1
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
	std::vector<int> values(start.back(), 1);  //����=start�����һ��Ԫ�أ�ȫ1����,�洢���е��ڽ���Ϣ  : start.back() = idx.size()

	const SE::SeCsr<int> m_neighbours_temp(start, idx, values); //const SE::SeCsr<int> m_neighbours(start,idx,values)
	const SE::SeCsr<int>* m_neighbours = &m_neighbours_temp; //ָ�������ֵ�ķ���: ������ͨ������ȡ��ַ��ֵ;����ʹ��new���帳ֵ
	//std::cout << m_neighbours->Size() << std::endl; //ָ�������->,��ͨ������.

	//2��PreparePreconditioner input
	//length = numverts
	//const SE::SeMatrix3f* diagonal = new SE::SeMatrix3f[numvert]{ SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };
	SE::SeMatrix3f* diagonal = new SE::SeMatrix3f[numvert]{ SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };

	//csrrange����start,�洢ÿ���ʵ��׸�offdiagnoal���б�csrrange[]�е�index  
	int* csrRanges = new int[numvert + 1]; //��ʾ��i�������ǰi-1������һ���ж����ھӣ�����˵Ӧ��=start
	int nn = 0;
	for (auto i = start.begin(); i != start.end(); i++) { csrRanges[nn] = (int)*i; nn++; }

	//std::cout << csrRanges[0] << std::endl;
	//�洢�ʵ�֮���������Ϣ,length = sum_{i}(neighbor[i]-1):���ж������ӵ��ھӸ������ܺ�
	//SE::SeMatrix3f* csrOffDiagonals = new SE::SeMatrix3f[start.back()]{ SE::SeMatrix3f::Identity(3),SE::SeMatrix3f::Identity(1),SE::SeMatrix3f::Identity(2),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity(),SE::SeMatrix3f::Identity() };
	SE::SeMatrix3f* csrOffDiagonals = new SE::SeMatrix3f[start.back()];
	
	//matrix3f��ֵ���㣺���ַ�ʽ�������Ҫ���ָ��Ԫ�أ���float list[9]��ֵ�����ֻ�����Խ��߸�ֵ��SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(a+1);
	//for (int a = 0; a < start.back(); a++)
	//{
	//	float r[9] = { 0,1,2,3,4,5,6,7,8 };
	//	SE::Float3 temp(1,2,3);
	//	SE::SeMatrix3f r1(r);//ע�������ǰ��и�ֵ�����ǰ���
	//	SE::SeMatrix3f r2 = SE::SeMatrix3f::Identity(1);
	//	SE::SeMatrix3f r3; 		r3.SetDiagonal(temp);
	//}

	for (int a = 0; a < start.back(); a++)
	{
		SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(0); //�ǶԽ�Ԫ��ָ����ͬ���½���Ĳ�ͬ
		csrOffDiagonals[a] = nows;
	}
	//�ο�����PSCC: �˴��ļ�����Ҫͨ����ײ�㷨���õ�,��ʾ��һʱ�̼�⵽��ײ����������
	const SE::EfSet* efSets = {}; //�Ȳ�����ef��ײ
	const SE::EeSet* eeSets = {};
	const SE::VfSet* vfSets = {};
	unsigned int* efCounts = new unsigned int[numedge + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 }; //�����б��ʼ��Ϊ0
	unsigned int* eeCounts = new unsigned int[numedge + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	unsigned int* vfCounts = new unsigned int[numvert + 1]{ 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

	//3��Preconditioning input
	//std::vector<SE::SeVec3fSimd> z;
	SE::SeVec3fSimd* z = new SE::SeVec3fSimd[numvert]; //��ʼΪ��,���ڴ洢�����Ľ��
	
	//residual �Ҳೣ����
	SE::SeVec3fSimd* residual = new SE::SeVec3fSimd[numvert]; //֮�����ֳɽ����ֵ��ʹ��const���壬������ɾ����
	for (int a = 0; a < numvert; a++)
	{
		SE::SeVec3fSimd nows(0.2, 0.1, 0.3, 0);
		residual[a] = nows;
	}
	int dim = 1;
	//4��prepare preconditioner 
	//SE::SeSchwarzPreconditioner* sh = nullptr; 
	SE::SeSchwarzPreconditioner sh;
	sh.m_positions = m_positions;
	sh.m_edges = m_edges;
	sh.m_faces = m_faces;
	sh.m_neighbours = m_neighbours; //��neighborΪ��ͨ����,�˴���ֵ���&

	sh.AllocatePrecoditioner(numvert, numedge, numface);
	std::cout << "AllocatePrecoditioner run success!" << std::endl;
	//
	///PCG����ʱ���ã�����PCGѭ���ÿ��ѭ���õ�ϵͳ����
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