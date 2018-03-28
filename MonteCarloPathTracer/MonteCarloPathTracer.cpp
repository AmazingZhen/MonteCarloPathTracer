#include "MonteCarloPathTracer.h"

#include <random>
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

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


glm::vec3 uniformSampleHemisphere(const glm::vec3 &normal) {
	glm::vec3 N, Nt, Nb;
	N = normal;

	if (std::fabs(N.x) > std::fabs(N.y))
		Nt = glm::vec3(N.z, 0, -N.x) / sqrtf(N.x * N.x + N.z * N.z);
	else
		Nt = glm::vec3(0, -N.z, N.y) / sqrtf(N.y * N.y + N.z * N.z);
	Nb = glm::cross(N, Nt);

	float r1 = distribution(generator);
	float r2 = distribution(generator);

	float sinTheta = sqrtf(1 - r1 * r1);
	float phi = 2 * M_PI * r2;
	float x = sinTheta * cosf(phi);
	float z = sinTheta * sinf(phi);
	
	glm::vec3 localRay(x, r1, z);

	return glm::vec3(
		localRay.x * Nb.x + localRay.y * N.x + localRay.z * Nt.x,
		localRay.x * Nb.y + localRay.y * N.y + localRay.z * Nt.y,
		localRay.x * Nb.z + localRay.y * N.z + localRay.z * Nt.z);
}

glm::vec3 MonteCarloPathTracer::trace(const Ray &ray, int depth)
{
	IntersectionInfo info;

	if (depth > options.maxDepth || !intersect(ray, info)) {
		return options.backgroundColor;
	}

	const Material &mater = *(info.mater);

	glm::vec3 surfaceLight = mater.ke;
	// Test ray casting.
	//static const glm::vec3 cols[3] = { { 0.6, 0.4, 0.1 },{ 0.1, 0.5, 0.3 },{ 0.1, 0.3, 0.7 } };
	//surfaceLight += info.localPoint.x * cols[0] + info.localPoint.y * cols[1] + info.localPoint.z * cols[2];
	//return surfaceLight;

	glm::vec3 directLight = directLighting(ray, info);
	surfaceLight += directLight;

	//float bias = 1e-4; // add some bias to the point from which we will be tracing
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
		glm::vec3 reflectLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1);
		
		float ior = mater.ni, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
		float cosi = -glm::dot(ray.d, normal);
		float k = 1 - eta * eta * (1 - cosi * cosi);
		glm::vec3 refractDir = ray.d * eta + normal * (eta *  cosi - sqrt(k));
		glm::normalize(refractDir);
		glm::vec3 refractLight = trace(Ray(info.point - normal * options.bias, refractDir), depth + 1);

		surfaceLight += reflectLight * fresneleffect
			+ refractLight * (1 - fresneleffect);
	}
	else if (glm::dot(mater.ks, mater.ks) > 0) {  // specular
		glm::vec3 reflectDir = ray.d - normal * 2.f * glm::dot(ray.d, normal);
		glm::normalize(reflectDir);
		glm::vec3 reflectLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1) * mater.ks;

		surfaceLight += reflectLight;
	}

	return surfaceLight;
}

bool MonteCarloPathTracer::intersect(const Ray & r, IntersectionInfo & info)
{
	return options.model ? options.model->intersect(r, info) : false;
}

glm::vec3 MonteCarloPathTracer::directLighting(const Ray &ray, IntersectionInfo &info)
{
	glm::vec3 directLight(0);

	if (glm::dot(info.mater->kd, info.mater->kd) <= 0.f) {
		return directLight;
	}

	for (uint32_t i = 0; i < options.lights.size(); ++i) {
		glm::vec3 lightDir, lightIntensity;

		options.lights[i].illuminate(info.point + info.normal * options.bias, lightDir, lightIntensity);

		Ray shadowRay(info.point + info.normal * options.bias, -lightDir);
		IntersectionInfo shadowInfo;
		intersect(shadowRay, shadowInfo);
		float shadowDistance = 0, lightDistance = 0;
		if (shadowInfo.hit) {
			shadowDistance = glm::dot(shadowInfo.point - shadowRay.o, shadowInfo.point - shadowRay.o);
			lightDistance = glm::dot(options.lights[i].pos - shadowRay.o, options.lights[i].pos - shadowRay.o);
		}

		if (!shadowInfo.hit || shadowDistance > lightDistance) {
			directLight += lightIntensity * std::max(0.f, glm::dot(info.normal, -lightDir)) * info.mater->kd;
		}
	}

	return directLight;
}

