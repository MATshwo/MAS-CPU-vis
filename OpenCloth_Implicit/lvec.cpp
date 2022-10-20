#include "lvec.h"
//������һЩ����Ĳ�����������ԭʼvector<glm::vec3>���������Ʋ����

float dot(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb) {
	float sum = 0;
	for (size_t i = 0; i < Va.v.size(); i++) {
		sum += glm::dot(Va.v[i], Vb.v[i]); //dot:���������һ������
	}
	return sum;
}
LargeVector<glm::vec3> operator*(const LargeVector<glm::mat3> other, const LargeVector<glm::vec3> v) {
	LargeVector<glm::vec3> tmp(v);
	for (size_t i = 0; i < v.v.size(); i++) {
		tmp.v[i] = other.v[i] * v.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator*(const LargeVector<glm::vec3> other, const LargeVector<glm::vec3> v) {
	LargeVector<glm::vec3> tmp(v);
	for (size_t i = 0; i < v.v.size(); i++) {
		tmp.v[i] = other.v[i] * v.v[i]; //���� = ����*����
	}
	return tmp;
}


LargeVector<glm::vec3> operator*(const float f, const LargeVector<glm::vec3> other) {
	LargeVector<glm::vec3> tmp(other);
	for (size_t i = 0; i < other.v.size(); i++) {
		tmp.v[i] = other.v[i] * f;
	}
	return tmp;
}
LargeVector<glm::mat3> operator*(const float f, const LargeVector<glm::mat3> other) {
	LargeVector<glm::mat3> tmp(other);
	for (size_t i = 0; i < other.v.size(); i++) {
		tmp.v[i] = other.v[i] * f;
	}
	return tmp;
}
LargeVector<glm::vec3> operator-(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb) {
	LargeVector<glm::vec3> tmp(Va);
	for (size_t i = 0; i < Va.v.size(); i++) {
		tmp.v[i] = Va.v[i] - Vb.v[i];
	}
	return tmp;
}
LargeVector<glm::mat3> operator-(const LargeVector<glm::mat3> Va, const LargeVector<glm::mat3> Vb) {
	LargeVector<glm::mat3> tmp(Va);
	for (size_t i = 0; i < Va.v.size(); i++) {
		tmp.v[i] = Va.v[i] - Vb.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator+(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb) {
	LargeVector<glm::vec3> tmp(Va);
	for (size_t i = 0; i < Va.v.size(); i++) {
		tmp.v[i] = Va.v[i] + Vb.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator/(const float f, const LargeVector<glm::vec3> v) {
	LargeVector<glm::vec3> tmp(v);
	for (size_t i = 0; i < v.v.size(); i++) {
		tmp.v[i] = v.v[i] / f;
	}
	return tmp;
}