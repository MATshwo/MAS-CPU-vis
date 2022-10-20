#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp> 
#include <glm/gtc/type_ptr.hpp>

template<typename T> class LargeVector {
private:
	std::vector<T> v;
public:
	LargeVector() {	}
	LargeVector(const LargeVector& other) {
		v.resize(other.v.size());
		memcpy(&v[0], &(other.v[0]), sizeof(other.v[0]) * other.v.size());
	}
	
	void resize(const int size){
		v.resize(size);
	}
	void clear(bool isIdentity = false) {
		memset(&v[0], 0, sizeof(T) * v.size());
		if (isIdentity) {
			for (unsigned int i = 0; i < v.size(); i++) {
				v[i] = T(1);
			}
		}
	}
	unsigned int size() {
		return v.size();
	}
	
	T &operator[](int index) {
		return v[index];
	}
	friend LargeVector<glm::vec3> operator*(const LargeVector<glm::mat3> other, const LargeVector<glm::vec3> f);
	friend LargeVector<glm::vec3> operator*(const float f, const LargeVector<glm::vec3> other);
	friend LargeVector<glm::vec3> operator-(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb);
	friend LargeVector<glm::vec3> operator*(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb);
	friend LargeVector<glm::vec3> operator+(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb);
	
	friend LargeVector<glm::mat3> operator*(const float f, const LargeVector<glm::mat3> other);
	friend LargeVector<glm::mat3> operator-(const LargeVector<glm::mat3> Va, const LargeVector<glm::mat3> Vb);
	friend LargeVector<glm::vec3> operator/(const float f, const LargeVector<glm::vec3> v);
	friend float dot(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb);
};

