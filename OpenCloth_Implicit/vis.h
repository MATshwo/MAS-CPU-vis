//��Ҫ���񣺵���MAS��implicit���ϣ��������µ�PCG��������

void SolveConjugateGradientMAS(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b,����ʩ�ߴ�Ԥ������sh);
void SolveConjugateGradientMAS(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b,����ʩ�ߴ�Ԥ������sh) {
	//Ԥ�����CG����� 
	//��ʼ��һ��z
	clock_t time_stt = clock();
	float i = 0;
	LargeVector<glm::vec3> r = (b - A * x);

	LargeVector<glm::vec3> z; //z0
	sh.precondition* (z, r, 1);//�����״�z
	LargeVector<glm::vec3> d = z; //p0


	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta = 0;
	float delta_old = 0;
	float delta_new = dot(r,z); //delta_new = dot(r,z); ��ʼ����
	float delta0 = delta_new;


	while (i<i_max && delta_new> EPS2 * delta0) {
		q = A * d;
		alpha = delta_new / dot(d, q);
		x = x + alpha * d;
		r = r - alpha * q;

		sh.precondition* (z, r, 1);//����µ�r����z
		delta_old = delta_new;
		delta_new = dot(z, r);

		beta = delta_new / delta_old; //���ϵ����
		d = r + beta * d; //d = z + beta * d; ����p
		i++;
	}
	//std::cout << "time use in IMPLICIT * 40 times is " << 1000 * (clock() - time_stt) / (double)CLOCKS_PER_SEC << "ms" << std::endl;
}
