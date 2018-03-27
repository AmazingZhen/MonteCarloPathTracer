#pragma once
#include <glm/gtc/matrix_transform.hpp>
#include <vector>

#include "Model.h"

const float M_PI = 3.1415926;
const float M_1_PI = 1.0 / M_PI;

struct PointLight
{
	glm::vec3 pos;
	glm::vec3 color;
	float intensity;

	PointLight(const glm::vec3 &p, const glm::vec3 &c = glm::vec3(1), const float &i = 10) :
		pos(p),
		color(c),
		intensity(i)
	{
	}
	// P: is the shaded point
	void illuminate(const glm::vec3 &P, glm::vec3 &lightDir, glm::vec3 &lightIntensity) const
	{
		lightDir = (P - pos);
		float r2 = glm::dot(lightDir, lightDir);
		// avoid division by 0
		lightIntensity = color * intensity / (4 * M_PI * r2);
		// printf("(%f, %f, %f)\n", lightIntensity.r, lightIntensity.g, lightIntensity.b);
	}
};

struct TracerOptions
{
	Model *model = 0;

	int width = 600;
	int height = 600;
	float fov = 70;
	float zNear = 0.1;
	float zFar = 1000.f;

	glm::vec3 eye = glm::vec3(0, 0, -2);
	glm::vec3 center = glm::vec3(0, 0, 0);
	glm::vec3 up = glm::vec3(0, 1, 0);

	int spp = 1;  // samples per pixel

	float bias = 0.0001;
	int maxDepth = 5;

	glm::vec3 backgroundColor = glm::vec3(0, 0, 0);

	std::vector<PointLight> lights;
};

class MonteCarloPathTracer
{
public:
	MonteCarloPathTracer(const TracerOptions &options);

	//void setViewport(int width, int height);
	void setCamera(const glm::vec3 &eye, const glm::vec3 &center, const glm::vec3 &up);
	void setPerspective(float fovy, float aspect, float zNear, float zFar);

	void render();
	void save();

private:
	Ray screenPointToRay(int r, int c);
	glm::vec3 trace(const Ray &ray, int depth);
	bool intersect(const Ray &r, IntersectionInfo &info);
	glm::vec3 directLighting(IntersectionInfo &info);


private:
	TracerOptions options;

	glm::mat4x4 camToWorld;
	glm::mat4x4 view;
	glm::mat4x4 perspective;

	std::vector<glm::vec3> color;
};
