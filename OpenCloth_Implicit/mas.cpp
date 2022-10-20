#include "mas.h"

namespace my
{
	
	//SE::SeSchwarzPreconditioner MAS(LargeVector<glm::mat3> sys_A, SE::SeSchwarzPreconditioner sh)
	//
	void MAS(LargeVector<glm::mat3> sys_A)
	{
		
		LARGE_INTEGER t0, tt, tc;
		QueryPerformanceFrequency(&tc);//tc��Ƶ�ʣ��������Ǽ���ϵͳ���Ĵ�������Ƶ��
		QueryPerformanceCounter(&t0);//��ʼʱ��
		int canPrintInfo = 0;
	
		

		//Loading data
		//Default Scale
		const int row_num = 301;
		const int col_num = 301;
		const int numvert = row_num * col_num; 
		const int numedge = (row_num-1)*(col_num)*3+(row_num)*2;
		const int numface = (row_num - 1) * (col_num) * 2;

		//vertex position
		SE::SeVec3fSimd* m_positions = new SE::SeVec3fSimd[numvert];
		int vcount = 0;
		for (int i = 0; i < row_num; i++)
		{
			for (int j = 0; j < col_num; j++)
			{
				SE::SeVec3fSimd temp(((float(i) / (row_num - 1)) * 2 - 1) * 2.0f, 4.0f + 1, ((float(j) / (col_num - 1)) * 4.0f)); 
				m_positions[vcount] = temp;
				vcount++;
			}
		}
		//Edge & Face
		SE::Int4* m_edges = new SE::Int4[numedge * 1.5]{ 0,0,0,0 };
		my::gen_edge(numvert, row_num, col_num, m_edges, m_positions);
		SE::Int4* m_faces = new SE::Int4[numface * 1.5];
		my::gen_face(numvert, row_num, col_num, m_faces, m_positions);
		//neighbor��Ϣ
		int start_temp[numvert + 1];
		int id_length = my::neibor_start(numvert, row_num, col_num, start_temp); // last value in start_temp 
		int* idx_temp = new int[id_length*1.5];
		my::neibor_idx(numvert, row_num, col_num, idx_temp);
		std::vector<int> start(start_temp, start_temp + numvert + 1);
		std::vector<int> idx(idx_temp, idx_temp + id_length);
		start.shrink_to_fit();
		idx.shrink_to_fit();
		std::vector<int> values(id_length, 1);
		SE::SeCsr<int> m_neighbours_temp(start, idx, values);
		SE::SeCsr<int>* m_neighbours = &m_neighbours_temp;

		//System matrix A value: diagonal & offdiagonal
		SE::SeMatrix3f* diagonal = new SE::SeMatrix3f[numvert * 1.5];
		SE::SeMatrix3f* csrOffDiagonals = new SE::SeMatrix3f[id_length * 1.5];
		int* csrRanges = new int[(numvert + 1) * 1.5];
		int nn = 0;


		OMP_PARALLEL
		{
			//��for������ò��л�
		#pragma omp for
			for (int a = 0; a < numvert; a++)
			{
			SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(1); //ֱ���õ�λ��
			diagonal[a] = nows;
			}
			
			for (auto i = start.begin(); i != start.end(); i++) { csrRanges[nn] = (int)*i; nn++; }
			for (int a = 0; a < id_length; a++)
			{
				SE::SeMatrix3f nows = SE::SeMatrix3f::Identity(0); //�ǶԽ�Ԫ��ָ����ͬ���½���Ĳ�ͬ
				csrOffDiagonals[a] = nows;
			}
		}
		
		////collision structure: �˴��޳���������ײ����
		////DCD���EF
		//SE::EfSet* efSets = {}; 
		////CCD���VF&EE
		//SE::EeSet* eeSets = {}; 
		//SE::VfSet* vfSets = {};
		//unsigned int* efCounts = new unsigned int[numedge + 1]{ 0 }; 
		//unsigned int* eeCounts = new unsigned int[numedge + 1]{ 0 };
		//unsigned int* vfCounts = new unsigned int[numvert + 1]{ 0 };


		//prepare&precompution preconditioner 
		sh.m_positions = m_positions;
		sh.m_edges = m_edges;
		sh.m_faces = m_faces;
		sh.m_neighbours = m_neighbours; 

		//input all parameters and initialize
		sh.AllocatePrecoditioner(numvert, numedge, numface); //���ø�һ��ʱ������
		//sh.PreparePreconditioner(diagonal, csrOffDiagonals, csrRanges,efSets,eeSets,vfSets,efCounts,eeCounts,vfCounts); //����ط�����ʱ��̫����
		sh.PreparePreconditioner(diagonal, csrOffDiagonals, csrRanges);
		
		delete m_positions;
		delete m_edges;
		delete m_faces;
		delete idx_temp;
		delete diagonal;
		delete csrRanges;
		delete csrOffDiagonals;
		//delete efCounts;
		//delete eeCounts;
		//delete vfCounts;

		if (canPrintInfo)//���Դ���bool������ȥ��ÿ��100֡�����
		{
			QueryPerformanceCounter(&tt);//��ǰʱ��
			//ͳ��AllocatePrecoditioner()&PreparePreconditioner()��ʱ
			std::cout << "MAS precomputation: " << (double)((tt.QuadPart - t0.QuadPart) * 1000.0 / tc.QuadPart) << "ms" << std::endl;
			//std::cout << "PreparePreconditioner run success!" << std::endl;
		}

		//return sh; //��main�ж���shΪȫ�ֱ������˴��Ͳ���Ҫ����ֵ���ٶ���30ms->3ms(˨Q)
	}


	//int Cmu(SE::SeSchwarzPreconditioner sh,LargeVector<glm::vec3>& z0,LargeVector<glm::vec3> sys_b)
	//{
	//	// execute matrix-vec multiplication
	//	//clock_t time_stt = clock();
	//	//3��Preconditioning input
	//	SE::SeVec3fSimd* z = new SE::SeVec3fSimd[sh.m_numVerts]; //���µĲ�
	//	SE::SeVec3fSimd* residual = new SE::SeVec3fSimd[sh.m_numVerts];
	//	for (int a = 0; a < sh.m_numVerts; a++)
	//	{
	//		SE::SeVec3fSimd nows(sys_b[a][0], sys_b[a][1], sys_b[a][2]); //��Ϊ����level����
	//		residual[a] = nows;
	//	}
	//	const int dim = 1;
	//	sh.Preconditioning(z, residual, dim);
	//	for (int k = 0; k < sh.m_numVerts; k++)
	//	{
	//		z0[k] = glm::vec3(z[k].x, z[k].y, z[k].z);
	//	}
	//	delete z;
	//	delete residual;
	//	return 0;
	//}
}