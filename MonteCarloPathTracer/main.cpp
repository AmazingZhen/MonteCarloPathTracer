#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2

#include "Model.h"
#include "MonteCarloPathTracer.h"

bool inter(const Ray &ray, const glm::vec3 &v0, const glm::vec3 &v1, const glm::vec3 &v2,
	float &t) {

	// compute plane's normal
	glm::vec3 v0v1 = v1 - v0;
	glm::vec3 v0v2 = v2 - v0;
	// no need to normalize
	glm::vec3 N = glm::cross(v0v1, v0v2); // N 
	float area2 = N.length();

	// Step 1: finding P

	// check if ray and plane are parallel ?
	float NdotRayDirection = glm::dot(N, ray.d);
	if (fabs(NdotRayDirection) < 1e-8) // almost 0 
		return false; // they are parallel so they don't intersect ! 

					  // compute d parameter using equation 2
	float d = glm::dot(N, v0);

	// compute t (equation 3)
	t = -(glm::dot(N, ray.o) + d) / NdotRayDirection;
	// check if the triangle is in behind the ray
	if (t < 0) return false; // the triangle is behind 

	// compute the intersection point using equation 1
	glm::vec3 P = ray.o + t * ray.d;

	// Step 2: inside-outside test
	glm::vec3 C; // vector perpendicular to triangle's plane 

				 // edge 0
	glm::vec3 edge0 = v1 - v0;
	glm::vec3 vp0 = P - v0;
	C = glm::cross(edge0, vp0);
	if (glm::dot(N, C) < 0) return false; // P is on the right side 

										  // edge 1
	glm::vec3 edge1 = v2 - v1;
	glm::vec3 vp1 = P - v1;
	C = glm::cross(edge1, vp1);
	if (glm::dot(N, C) < 0)  return false; // P is on the right side 

										   // edge 2
	glm::vec3 edge2 = v0 - v2;
	glm::vec3 vp2 = P - v2;
	C = glm::cross(edge2, vp2);
	if (glm::dot(N, C) < 0) return false; // P is on the right side; 

	return true; // this ray hits the triangle 
}

void testTriangleIntersect() {
	glm::vec3 v0(0, 0, 0), v1(1, 0, 0), v2(0, 1, 0);
	glm::vec3 o(0, 0, 1), d(0, 0, 1);

	float t = 0.f;

	while (std::cin >> o.x >> o.y >> o.z) {
		Ray ray(o, d);

		bool res = inter(ray, v0, v1, v2, t);
		std::cout << res << std::endl;
		if (res) {
			glm::vec3 point = ray.o + t * ray.d;
			printf("%f :(%f, %f, %f)\n", t, point.x, point.y, point.z);
		}
	}
}

int main(int argc, char *argv[]) {
	//testTriangleIntersect();

	Model m;

	m.loadModel("../scenes/scene01.obj");

	TracerOptions options;
	options.model = &m;

	// camera 
	glm::vec3 center = (m.getBoundingBox().bounds[0] + m.getBoundingBox().bounds[1]) / 2.f;
	float scale = glm::length(m.getBoundingBox().bounds[0] - m.getBoundingBox().bounds[1]) / 2;
	options.eye = center + glm::vec3(0, -.2f * scale, 1.5 * scale);
	options.center = center + glm::vec3(0, -.2f * scale, 0);

	PointLight light(center);
	//PointLight light(glm::vec3(-1, 9.8, 1));
	options.lights.push_back(light);

	MonteCarloPathTracer tracer(options);

	tracer.render();
	tracer.save();

	return 0;
}
