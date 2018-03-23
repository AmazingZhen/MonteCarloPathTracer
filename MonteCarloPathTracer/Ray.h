#pragma once

#include <glm/common.hpp>

#include "Mesh.h"

struct Ray {
	Ray(const glm::vec3 &o_, const glm::vec3 &d_) {
		o = o_;
		d = d_;

		invd = 1.f / d;
		sign[0] = (invd.x < 0);
		sign[1] = (invd.y < 0);
		sign[2] = (invd.z < 0);
	}

	glm::vec3 o;
	glm::vec3 d;
	glm::vec3 invd;
	int sign[3];
};