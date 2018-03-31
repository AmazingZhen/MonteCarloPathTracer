#include "MonteCarloPathTracer.h"

#include <random>
std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0, 1);

inline float clamp(float x) { return x<0 ? 0 : x>1 ? 1 : x; }

bool russianRoulette(float probability, float& survival)
{
	//srand(0);
	survival = distribution(generator);
	if (survival > probability) return true;
	return false;
}

MonteCarloPathTracer::MonteCarloPathTracer(const TracerOptions &options_)
{
	this->options = options_;

	setCamera(options.eye, options.center, options.up);
	setPerspective(options.fov, (float)options.width / options.height, options.zNear, options.zFar);

	color = std::vector<glm::vec3>(options.width * options.height, options.backgroundColor);
}

void MonteCarloPathTracer::setCamera(const glm::vec3 &eye, const glm::vec3 &center, const glm::vec3 &up_)
{
	view = glm::lookAt(eye, center, up_);
	camToWorld = glm::inverse(view);

	glm::vec3 direction = glm::normalize(center - eye);
	glm::vec3 right = glm::normalize(glm::cross(direction, up_));
	glm::vec3 up = glm::normalize(glm::cross(right, direction));

	float aspect = (float)options.height / (float)options.width;

	view_x = right * 2.f * tan(options.fov * M_PI / 360.f);
	view_y = up * 2.f * tan(options.fov * aspect * M_PI / 360.f);
	view_z = direction;
}

void MonteCarloPathTracer::setPerspective(float fovy, float aspect, float zNear, float zFar)
{
	perspective = glm::perspective(fovy, aspect, zNear, zFar);
}

void MonteCarloPathTracer::render()
{
	for (int iteration = 1; iteration <= options.spp; iteration++) {

#pragma omp parallel for schedule(dynamic, 1)
		for (int row = 0; row < options.height; row++) {
			fprintf(stderr, "\rRendering (%d iterations) %5.2f%%", iteration, 100.* row / (options.height - 1));

			for (int col = 0; col < options.width; col++) {
				int i = (options.height - row - 1) * options.width + col;
				//int i = row * options.width + col;
				Ray ray = screenPointToRay(row, col);
				glm::vec3 radiance = trace(ray, 0);

				color[i] = (color[i] * (iteration - 1.f) + radiance) / (float)iteration;
			}
		}

		save(iteration);
	}
}

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

void MonteCarloPathTracer::save(int iteration)
{
	char file[256];

	if (options.iterationsNeedSave.find(iteration) != options.iterationsNeedSave.end()) {
		sprintf_s(file, 256, "%s/image_iterations_%d.ppm", options.saveDir, iteration);
	}
	else if (iteration % 100 == 0) {
		sprintf_s(file, 256, "%s/image.ppm", options.saveDir);
	}

	FILE *f = fopen(file, "w");         // Write image to PPM file.
	if (f) {
		fprintf(f, "P3\n%d %d\n%d\n", options.width, options.height, 255);
		for (int i = 0; i < options.width * options.height; i++)
			fprintf(f, "%d %d %d ", toInt(clamp(color[i].x)), toInt(clamp(color[i].y)), toInt(clamp(color[i].z)));
		fclose(f);
	}
}

Ray MonteCarloPathTracer::screenPointToRay2(int row, int col)
{
	/*
		Camera coordinates.
	*/

	float x = col + (1 - 2 * distribution(generator));
	float y = row + (1 - 2 * distribution(generator));

	float imageAspectRatio = (float)options.width / options.height; // assuming width > height 
	float Px = (2 * ((x + 0.5) / options.width) - 1) * std::tan(options.fov / 2 * M_PI / 180)* imageAspectRatio;
	float Py = (1 - 3 * ((y + 0.5) / options.height) * std::tan(options.fov / 2 * M_PI / 180));
	
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


glm::vec3 uniformSampleHemisphere2(const glm::vec3 &normal) {
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

	return glm::normalize(glm::vec3(
		localRay.x * Nb.x + localRay.y * N.x + localRay.z * Nt.x,
		localRay.x * Nb.y + localRay.y * N.y + localRay.z * Nt.y,
		localRay.x * Nb.z + localRay.y * N.z + localRay.z * Nt.z));
}

glm::vec3 uniformSampleHemisphere(const glm::vec3 &normal) {
	float phi, theta;
	float r1 = distribution(generator), r2 = distribution(generator);

	phi = r1 * 2 * M_PI;
	theta = asinf(sqrtf(r2));
	glm::vec3 sample(sinf(theta) * cosf(phi), cosf(theta), sinf(theta) * sinf(phi));

	glm::vec3 front, right, up = normal;

	if (fabs(up.x) > fabs(up.y))
		front = glm::vec3(up.z, 0, -up.x);
	else
		front = glm::vec3(0, -up.z, up.y);

	front = normalize(front);
	right = cross(up, front);

	return glm::normalize(sample.x * right + sample.y * up + sample.z * front);
}

glm::vec3 importanceSample(const glm::vec3 &direction, float n) {
	float phi, theta;
	float r1 = distribution(generator), r2 = distribution(generator);

	phi = r1 * 2 * M_PI;
	theta = acosf(powf(r2, 1 / (n + 1)));
	glm::vec3 sample(sinf(theta) * cosf(phi), cosf(theta), sinf(theta) * sinf(phi));

	glm::vec3 front, right, up = direction;

	if (fabs(up.x) > fabs(up.y))
		front = glm::vec3(up.z, 0, -up.x);
	else
		front = glm::vec3(0, -up.z, up.y);

	front = normalize(front);
	right = cross(up, front);

	return glm::normalize(sample.x * right + sample.y * up + sample.z * front);
}

Ray MonteCarloPathTracer::screenPointToRay(int row, int col)
{
	float x = (col + (1 - 2 * distribution(generator))) / options.width;
	float y = (row + (1 - 2 * distribution(generator))) / options.height;

	glm::vec3 direction = view_z + (x - 0.5f) * view_x
		+ (y - 0.5f) * view_y;

	return Ray(options.eye, direction);
}

glm::vec3 MonteCarloPathTracer::trace(const Ray &ray, int depth)
{
	IntersectionInfo info;

	if (!intersect(ray, info)) {
		return options.backgroundColor;
	}

	const Material &mater = *(info.mater);

	glm::vec3 surfaceLight = mater.ke;

	if (depth > options.maxDepth) {
		return surfaceLight;
	}

	// Test ray casting.
	//static const glm::vec3 cols[3] = { { 0.6, 0.4, 0.1 },{ 0.1, 0.5, 0.3 },{ 0.1, 0.3, 0.7 } };
	//surfaceLight += info.localPoint.x * cols[0] + info.localPoint.y * cols[1] + info.localPoint.z * cols[2];
	//return surfaceLight;

	//glm::vec3 directLight = directLighting(ray, info);
	//surfaceLight += directLight;

	bool inside = false;
	glm::vec3 normal = info.normal;
	if (glm::dot(ray.d, info.normal) > 0) {
		normal = -normal;
		inside = true;
	}

	if (mater.ni != 1.0) {  // we think this is glass
		glm::vec3 reflectDir = ray.d - normal * 2.f * glm::dot(ray.d, normal);
		glm::normalize(reflectDir);
		
		float ior = mater.ni, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
		float cosi = -glm::dot(ray.d, normal);
		float k = 1 - eta * eta * (1 - cosi * cosi);

		float facingratio = -glm::dot(ray.d, normal);

		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);

		glm::vec3 refractDir = ray.d * eta + normal * (eta *  cosi - sqrt(k));
		glm::normalize(refractDir);

		glm::vec3 reflectLight(0);
		glm::vec3 refractLight(0);

		if (depth > 2) {
			float p = 0.5 * fresneleffect, survival;

			if (russianRoulette(p, survival)) {
				refractLight = trace(Ray(info.point - normal * options.bias, refractDir), depth + 1) * (1 - fresneleffect / 1 - p);
			}
			else {
				reflectLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1) * (fresneleffect / p);
			}
		}
		else {
			reflectLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1) * fresneleffect;
			refractLight = trace(Ray(info.point - normal * options.bias, refractDir), depth + 1) * (1 - fresneleffect);
		}

		surfaceLight += reflectLight + refractLight;
	}
	else {
		float spec = glm::dot(mater.ks, glm::vec3(1));
		float total = spec + glm::dot(mater.kd, glm::vec3(1));
		float p = spec / total, survival;

		if (russianRoulette(p, survival)) {
			glm::vec3 reflectDir = uniformSampleHemisphere(normal);
			glm::vec3 diffuseLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1) * mater.kd * (1 / 1 - p);

			surfaceLight += diffuseLight;
		}
		else  {  // specular
			glm::vec3 perfectReflectDir = glm::normalize(ray.d - normal * 2.f * glm::dot(ray.d, normal));
			glm::vec3 reflectDir = importanceSample(perfectReflectDir, mater.ns);

			glm::vec3 reflectLight = trace(Ray(info.point + normal * options.bias, reflectDir), depth + 1) * mater.ks * (1 / p);

			surfaceLight += reflectLight;
		}
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

