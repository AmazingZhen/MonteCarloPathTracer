#pragma once

#include <vector>
#include "Mesh.h"

struct KDNode
{
	KDNode(){};

	static KDNode* build(std::vector<Triangle*> &tris, int depth) {
		if (tris.empty()) {
			return 0;
		}
		
		KDNode *node = new KDNode();

		if (tris.size() == 1) {
			node->tri = tris[0];
			node->box = tris[0]->box;
			return node;
		}

		// Find longest axis;
		node->box = tris[0]->box;
		for (int i = 1; i < tris.size(); i++) {
			node->box.expand(tris[i]->box);
		}
		int axis = node->box.getLongestAxis();
		node->axis = axis;

		// Find mid value
		std::sort(tris.begin(), tris.end(), [axis](Triangle* a, Triangle* b) {
			return a->box.getMidPoint()[axis] < b->box.getMidPoint()[axis];
		});
		int midPos = ((tris.size() % 2) ? tris.size() / 2 : tris.size() / 2 - 1);
		node->tri = tris[midPos];

		std::vector<Triangle*> left_tris(tris.begin(), tris.begin() + midPos),
			right_tris(tris.begin() + midPos + 1, tris.end());

		//printf("level %d : %d + %d = %d\n", depth,left_tris.size(), right_tris.size(), tris.size());

		if (!left_tris.empty()) {
			node->left = build(left_tris, depth + 1);
		}
		
		if (!right_tris.empty()) {
			node->right = build(right_tris, depth + 1);
		}

		return node;
	}

	static bool intersect(KDNode *node, const Ray &ray, float &t_min, IntersectionInfo &info) {
		if (!node || !node->box.intersect(ray)) {
			return false;
		}
		
		bool intersectTriangle = false;

		float t, u, v;
		if (node->tri->intersect(ray, t, u, v) && t < t_min) {
			t_min = t;
			node->tri->getSurfaceProperties(ray, t, u, v, info);

			intersectTriangle = true;
		}

		intersectTriangle |= intersect(node->left, ray, t_min, info);
		intersectTriangle |= intersect(node->right, ray, t_min, info);

		return intersectTriangle;
	}

	int axis = 0;

	AABB box;
	Triangle* tri = 0;

	KDNode *left = 0;
	KDNode *right = 0;
};
