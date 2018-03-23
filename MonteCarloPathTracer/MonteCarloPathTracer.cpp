#include "MonteCarloPathTracer.h"


inline float clamp(float x) { return x<0 ? 0 : x>1 ? 1 : x; }

MonteCarloPathTracer::MonteCarloPathTracer(const TracerOptions &options_)
{
	this->options = options_;

	setCamera(options.eye, options.center, options.up);
	setPerspective(options.fov, (float)options.width / options.height, options.zNear, options.zFar);

	color = std::vector<glm::vec3>(options.width * options.height, options.backgroundColor);
}

void MonteCarloPathTracer::setCamera(const glm::vec3 &eye, const glm::vec3 &center, const glm::vec3 &up)
{
	view = glm::lookAt(eye, center, up);
	camToWorld = glm::inverse(view);
}

void MonteCarloPathTracer::setPerspective(float fovy, float aspect, float zNear, float zFar)
{
	perspective = glm::perspective(fovy, aspect, zNear, zFar);
}

void MonteCarloPathTracer::render()
{
//#pragma omp parallel for schedule(dynamic, 1)
	for (int row = 0; row < options.height; row++) {
		fprintf(stderr, "\rRendering (%d spp) %5.2f%%", options.spp, 100.* row / (options.height - 1));
		
		for (int col = 0; col < options.width; col++) {
			int i = (options.height - row - 1) * options.width + col;
			glm::vec3 radiance = glm::vec3(0, 0, 0);

			for (int s = 0; s < options.spp; s++) {
				Ray ray = screenPointToRay(row, col);
				/*
				printf("(%f, %f, %f), (%f, %f, %f)\n",
					ray.o.x, ray.o.y, ray.o.z,
					ray.d.x, ray.d.y, ray.d.z);
				int t;
				std::cin >> t;
				*/

				radiance += traceRadiance(ray, 0);
			}
			radiance /= options.spp;
			color[i] += glm::vec3(clamp(radiance.x), clamp(radiance.y), clamp(radiance.z));
		}
	}
}

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

void MonteCarloPathTracer::save()
{
	FILE *f = fopen("image.ppm", "w");         // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", options.width, options.height, 255);
	for (int i = 0; i< options.width * options.height; i++)
		fprintf(f, "%d %d %d ", toInt(color[i].x), toInt(color[i].y), toInt(color[i].z));
}

Ray MonteCarloPathTracer::screenPointToRay(int row, int col)
{
	/*
		Camera coordinates.
	*/

	float imageAspectRatio = (float)options.width / options.height; // assuming width > height 
	float Px = (2 * ((col + 0.5) / options.width) - 1) * std::tan(options.fov / 2 * M_PI / 180) * imageAspectRatio;
	float Py = (1 - 2 * ((row + 0.5) / options.height) * std::tan(options.fov / 2 * M_PI / 180));
	
	glm::vec3 rayOrigin(0);
	glm::vec3 rayDirection = glm::vec3(Px, Py, -1) - rayOrigin; // note that this just equal to Vec3f(Px, Py, -1); 
	rayDirection = normalize(rayDirection); // it's a direction so don't forget to normalize

	/*
		World coordinates.
	*/

	return Ray(glm::vec3(camToWorld * glm::vec4(rayOrigin, 1))
		, glm::vec3(camToWorld * glm::vec4(rayDirection, 0)));
}

glm::vec3 MonteCarloPathTracer::traceRadiance(const Ray &ray, int depth)
{
	//glm::vec3 hitColor = (r.d + glm::vec3(1)) * 0.5f;
	//return hitColor;

	IntersectionInfo info;

	if (!intersect(ray, info)) {
		return options.backgroundColor;
	}

	glm::vec3 directLight(0);
	for (uint32_t i = 0; i < options.lights.size(); ++i) {
		glm::vec3 lightDir, lightIntensity;

		options.lights[i].illuminate(info.point, lightDir, lightIntensity);

		Ray shadowRay(info.point, lightDir);
		IntersectionInfo shadowInfo;
		if (!intersect(shadowRay, shadowInfo)) {
			//printf("lightIntensity : (%f, %f, %f)\n", lightIntensity.r, lightIntensity.g, lightIntensity.b);
			//directLight += lightIntensity * std::max(0.f, glm::dot(info.normal, -lightDir)) * info.mater->kd;
			directLight += lightIntensity * info.mater->kd;

			//printf("directLight : (%f, %f, %f)\n", directLight.r, directLight.g, directLight.b);

			//int t;
			//std::cin >> t;
		}
	}

	return directLight;
}

bool MonteCarloPathTracer::intersect(const Ray & r, IntersectionInfo & info)
{
	return options.model ? options.model->intersect(r, info) : false;
}

