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
			//int i = (options.height - row - 1) * options.width + col;
			int i = row * options.width + col;
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

				radiance += trace(ray, 0);
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

float mix(const float &a, const float &b, const float &mix)
{
	return b * mix + a * (1 - mix);
}

glm::vec3 MonteCarloPathTracer::trace(const Ray &ray, int depth)
{
	IntersectionInfo info;

	if (!intersect(ray, info)) {
		return options.backgroundColor;
	}

	glm::vec3 surfaceLight(0);
	// Test ray casting.
	//static const glm::vec3 cols[3] = { { 0.6, 0.4, 0.1 },{ 0.1, 0.5, 0.3 },{ 0.1, 0.3, 0.7 } };
	//surfaceLight += info.localPoint.x * cols[0] + info.localPoint.y * cols[1] + info.localPoint.z * cols[2];
	//return surfaceLight;

	glm::vec3 directLight = directLighting(info);
	surfaceLight += directLight;
	if (depth == options.maxDepth) {
		return surfaceLight;
	}

	const Material &mater = *(info.mater);

	float bias = 1e-4; // add some bias to the point from which we will be tracing
	bool inside = false;
	glm::vec3 normal = info.normal;
	if (glm::dot(ray.d, info.normal) > 0) {
		normal = -normal;
		inside = true;
	}

	if (mater.ni != 1.0) {  // we think this is glass
		float facingratio = -glm::dot(ray.d, normal);

		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
		
		glm::vec3 reflectDir = ray.d - normal * 2.f * glm::dot(ray.d, normal);
		glm::normalize(reflectDir);
		glm::vec3 reflectLight = trace(Ray(info.point + normal * bias, reflectDir), depth + 1);
		
		float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
		float cosi = -glm::dot(ray.d, normal);
		float k = 1 - eta * eta * (1 - cosi * cosi);
		glm::vec3 refractDir = ray.d * eta + normal * (eta *  cosi - sqrt(k));
		glm::normalize(refractDir);
		glm::vec3 refractLight = trace(Ray(info.point - normal * bias, refractDir), depth + 1);

		surfaceLight += reflectLight * fresneleffect
			+ refractLight * (1 - fresneleffect);
	}
	else {
		glm::vec3 reflectDir = ray.d - normal * 2.f * glm::dot(ray.d, normal);
		glm::normalize(reflectDir);
		glm::vec3 reflectLight = trace(Ray(info.point + normal * bias, reflectDir), depth + 1) * mater.ks;

		surfaceLight += reflectLight;
	}

	return surfaceLight + mater.ke;
}

bool MonteCarloPathTracer::intersect(const Ray & r, IntersectionInfo & info)
{
	return options.model ? options.model->intersect(r, info) : false;
}

glm::vec3 MonteCarloPathTracer::directLighting(IntersectionInfo & info)
{
	glm::vec3 directLight(0);

	for (uint32_t i = 0; i < options.lights.size(); ++i) {
		glm::vec3 lightDir, lightIntensity;

		options.lights[i].illuminate(info.point, lightDir, lightIntensity);

		Ray shadowRay(info.point, lightDir);
		IntersectionInfo shadowInfo;
		if (!intersect(shadowRay, shadowInfo)) {
			//printf("lightIntensity : (%f, %f, %f)\n", lightIntensity.r, lightIntensity.g, lightIntensity.b);

			//directLight += info.localPoint.x * cols[0] + info.localPoint.y * cols[1] + info.localPoint.z * cols[2];
			directLight += lightIntensity * std::max(0.f, glm::dot(info.normal, -lightDir)) * info.mater->kd;
			//directLight += lightIntensity * info.mater->kd;

			//printf("directLight : (%f, %f, %f)\n", directLight.r, directLight.g, directLight.b);
		}
	}

	return directLight;
}

