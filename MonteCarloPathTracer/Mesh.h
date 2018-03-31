#pragma once

#include <glm/glm.hpp>
#include <iostream>

#include "AABB.h"

struct Material
{
	std::string name;

	glm::vec3 ka;  // ambient
	glm::vec3 kd;  // diffuse
	glm::vec3 ks;  // specular
	glm::vec3 ke;  // emission

	float ns;  // shininess
	float ni;  // Index Of Refraction

	glm::vec3 tf;  // Transparency

	void print() {
		printf("name:%s\n", name.c_str());
		printf("ka: %f, %f, %f\n", ka.x, ka.y, ka.z);
		printf("kd: %f, %f, %f\n", kd.x, kd.y, kd.z);
		printf("ks: %f, %f, %f\n", ks.x, ks.y, ks.z);
		printf("ke: %f, %f, %f\n", ke.x, ke.y, ke.z);

		printf("ns: %f\n", ns);
		printf("ni: %f\n", ni);

		printf("tf: %f, %f, %f\n", tf.x, tf.y, tf.z);
		printf("\n");
	}
};

struct Mesh {
	std::vector<glm::vec3> vertices;
	std::vector<unsigned int> indices;
	std::vector<glm::vec3> normals;

	// unsigned int numIndices;
	unsigned int materialIndex;
};


struct IntersectionInfo {
	bool hit = false;

	glm::vec3 point = glm::vec3(0);
	glm::vec3 normal = glm::vec3(0);

	glm::vec3 localPoint = glm::vec3(0);

	Material* mater = 0;
};

struct Triangle {
	Triangle(unsigned int indices_[3], Mesh *mesh_, Material *mater_) {
		indices[0] = indices_[0];
		indices[1] = indices_[1];
		indices[2] = indices_[2];

		mesh = mesh_;
		mater = mater_;

		glm::vec3 low, high;
		glm::vec3& pt0 = mesh->vertices[indices[0]],
			pt1 = mesh->vertices[indices[1]],
			pt2 = mesh->vertices[indices[2]];

		low.x = std::min(pt0.x, std::min(pt1.x, pt2.x));
		low.y = std::min(pt0.y, std::min(pt1.y, pt2.y));
		low.z = std::min(pt0.z, std::min(pt1.z, pt2.z));

		high.x = std::max(pt0.x, std::max(pt1.x, pt2.x));
		high.y = std::max(pt0.y, std::max(pt1.y, pt2.y));
		high.z = std::max(pt0.z, std::max(pt1.z, pt2.z));

		box = AABB(low, high);
	}

	bool intersect(const Ray &ray, float &t, float &u, float &v) {
		const glm::vec3 &v0 = mesh->vertices[indices[0]];
		const glm::vec3 &v1 = mesh->vertices[indices[1]];
		const glm::vec3 &v2 = mesh->vertices[indices[2]];

		glm::vec3 edge1, edge2, h, s, q;
		float a, f;
		edge1 = v1 - v0;
		edge2 = v2 - v0;
		h = glm::cross(ray.d, edge2);
		a = glm::dot(edge1, h);
		if (a > -1e-5 && a < 1e-5)
			return false;
		
		f = 1 / a;
		s = ray.o - v0;
		u = f * (glm::dot(s, h));
		if (u < 0.0 || u > 1.0)
			return false;
		q = glm::cross(s, edge1);
		v = f * glm::dot(ray.d,q);
		if (v < 0.0 || u + v > 1.0)
			return false;
		// At this stage we can compute t to find out where the intersection point is on the line.
		t = f * glm::dot(edge2, q);
		if (t > 1e-5) // ray intersection
		{
			return true;
		}
		else // This means that there is a line intersection but not a ray intersection.
			return false;
	}

	void getSurfaceProperties(
		const Ray &ray,
		float t,
		float u,
		float v,
		IntersectionInfo &info) const
	{
		info.hit = true;

		info.point = ray.o + t * ray.d;

		/*
		const glm::vec3 &v0 = mesh->vertices[indices[0]];
		const glm::vec3 &v1 = mesh->vertices[indices[1]];
		const glm::vec3 &v2 = mesh->vertices[indices[2]];

		printf("v0 : (%f, %f, %f)\n", v0.x, v0.y, v0.z);
		printf("v1 : (%f, %f, %f)\n", v1.x, v1.y, v1.z);
		printf("v2 : (%f, %f, %f)\n", v2.x, v2.y, v2.z);
		printf("point : (%f, %f, %f)\n", info.point.x, info.point.y, info.point.z);
		*/

		info.localPoint.x = (1 - u - v);
		info.localPoint.y = u;
		info.localPoint.z = v;

		const std::vector<glm::vec3> &N = mesh->normals;

		// vertex normal
		const glm::vec3 &n0 = N[indices[0]];
		const glm::vec3 &n1 = N[indices[1]];
		const glm::vec3 &n2 = N[indices[2]];
		info.normal = glm::normalize((1 - u - v) * n0 + u * n1 + v * n2);

		info.mater = mater;
	}

	unsigned int indices[3];

	Mesh *mesh;
	Material *mater;

	AABB box;
};
