#pragma once

#include <glm/common.hpp>
#include <algorithm>

#include "Ray.h"

// Axis-Align Bounding Box
struct AABB
{
public:
	AABB() { bounds[0] = bounds[1] = glm::vec3(0); };
	AABB(glm::vec3 low, glm::vec3 high) {
		bounds[0] = low;
		bounds[1] = high;
	}

	void expand(const AABB &box) {
		glm::vec3 low, high;

		low.x = std::min(bounds[0].x, box.bounds[0].x);
		low.y = std::min(bounds[0].y, box.bounds[0].y);
		low.z = std::min(bounds[0].z, box.bounds[0].z);

		high.x = std::max(bounds[1].x, box.bounds[1].x);
		high.y = std::max(bounds[1].y, box.bounds[1].y);
		high.z = std::max(bounds[1].z, box.bounds[1].z);

		bounds[0] = low;
		bounds[1] = high;
	};

	glm::vec3 getMidPoint() {
		return (bounds[0] + bounds[1]) / 2.f;
	}

	int getLongestAxis() {
		glm::vec3 diff = bounds[1] - bounds[0];
		float longestAxisValue = std::max(diff.x,
			std::max(diff.y, diff.z));

		for (int i = 0; i < diff.length(); i++) if (diff[i] == longestAxisValue) {
			return i;
		}
	}

	bool intersect(const Ray &r) const
	{
		float tmin, tmax, tymin, tymax, tzmin, tzmax;

		tmin = (bounds[r.sign[0]].x - r.o.x) * r.invd.x;
		tmax = (bounds[1 - r.sign[0]].x - r.o.x) * r.invd.x;
		tymin = (bounds[r.sign[1]].y - r.o.y) * r.invd.y;
		tymax = (bounds[1 - r.sign[1]].y - r.o.y) * r.invd.y;

		if ((tmin > tymax) || (tymin > tmax))
			return false;
		if (tymin > tmin)
			tmin = tymin;
		if (tymax < tmax)
			tmax = tymax;

		tzmin = (bounds[r.sign[2]].z - r.o.z) * r.invd.z;
		tzmax = (bounds[1 - r.sign[2]].z - r.o.z) * r.invd.z;

		if ((tmin > tzmax) || (tzmin > tmax))
			return false;
		if (tzmin > tmin)
			tmin = tzmin;
		if (tzmax < tmax)
			tmax = tzmax;

		return true;
	}

	glm::vec3 bounds[2];
};
