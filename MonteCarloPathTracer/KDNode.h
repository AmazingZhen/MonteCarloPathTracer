#pragma once

#include <vector>
#include "Mesh.h"

struct KDNode
{
	KDNode(){};

	static KDNode* build(const std::vector<Triangle*>& tris, int depth) {
		if (tris.empty()) {
			return 0;
		}
		
		KDNode *node = new KDNode();
		node->tris = tris;
		node->box = tris[0]->box;

		if (tris.size() == 1) {
			return node;
		}

		for (int i = 1; i < tris.size(); i++) {
			node->box.expand(tris[i]->box);
		}

		glm::vec3 midPoint(0);
		for (int i = 0; i < tris.size(); i++) {
			midPoint += tris[i]->box.getMidPoint();
		}
		midPoint /= tris.size();

		std::vector<Triangle*> left_tris, right_tris;
		int axis = node->box.getLongestAxis();
		for (Triangle* tri : tris) {
			if (midPoint[axis] >= tri->box.getMidPoint()[axis]) {
				left_tris.push_back(tri);
			}
			else {
				right_tris.push_back(tri);
			}
		}

		if (left_tris.size() == 0 && right_tris.size() > 0) {
			left_tris = right_tris;
		}

		if (right_tris.size() == 0 && left_tris.size() > 0) {
			right_tris = left_tris;
		}

		float matches = 0;
		for (Triangle* tl : left_tris) {
			for (Triangle* tr : right_tris) if (tl == tr){
				matches++;
			}
		}

		if (matches / left_tris.size() < 0.5 && matches / right_tris.size() < 0.5) {
			node->left = build(left_tris, depth + 1);
			node->right = build(right_tris, depth + 1);
		}

		return node;
	}

	static bool intersect(KDNode *node, const Ray &ray, float &t, float &tmin, IntersectionInfo &info) {
		if (!node) {
			return false;
		}
		
		if (node->box.intersect(ray)) {
			bool hitTriangle = false;
			glm::vec3 intersectPoint, intersectNormal, localIntersectPoint;

			if (node->left || node->right) {
				bool intersectLeft = intersect(node->left, ray, t, tmin, info);
				bool intersectRight = intersect(node->right, ray, t, tmin, info);

				return intersectLeft || intersectRight;
			}
			else {  // reach a leaf node
				float u, v;

				for (Triangle *tri : node->tris) {
					if (tri->intersect(ray, t, u, v)) {
						hitTriangle = true;
						info.hit = true;
						tmin = t;
						intersectPoint = info.point;
						intersectNormal = info.normal;
						localIntersectPoint = info.localPoint;
					}
				}

				if (hitTriangle) {
					info.hit = true;

					info.point = intersectPoint;
					info.normal = intersectNormal;
					info.localPoint = localIntersectPoint;

					return true;
				}

				return false;
			}
		}

		return false;
	}

	AABB box;
	std::vector<Triangle*> tris;

	KDNode *left = 0;
	KDNode *right = 0;
};
