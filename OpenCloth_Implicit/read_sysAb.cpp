#include "read_sysAb.h"

namespace my
{ 
	float* read_data(char* path, int length)
	{
		float *list = new float[length];
		std::ifstream in(path);
		for (int i = 0; i < length; i++) { in >> list[i]; }
		//for (int i = 0; i < length; i++) { std::cout << list[i] << std::endl; }
		in.close();
		return list;
	}
}



