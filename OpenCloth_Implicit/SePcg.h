#pragma once


void SolveConjugateGradient(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b);

void SolveConjugateGradient(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b)
{
	float i = 0;
	LargeVector<glm::vec3> r = b - A * x;
	LargeVector<glm::vec3> d = r;
	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta = 0;
	float delta_old = 0;
	float delta_new = dot(r, r);
	float delta0 = delta_new;
	while (i<i_max && delta_new> EPS2) {
		q = A * d;
		alpha = delta_new / dot(d, q);
		x = x + alpha * d;
		r = r - alpha * q;
		delta_old = delta_new;
		delta_new = dot(r, r);
		beta = delta_new / delta_old;
		d = r + beta * d;
		i++;
	}
	}
}