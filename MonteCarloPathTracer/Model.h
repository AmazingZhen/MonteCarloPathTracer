#pragma once

#include <iostream>
#include <vector>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <assimp/Importer.hpp>
#include <glm/common.hpp>

#include "Mesh.h"
#include "KDNode.h"

class Model
{
public:
	Model() {};

	~Model() {};

	bool loadModel(const std::string& filename);
	bool intersect(const Ray & r, IntersectionInfo & info);

	AABB getBoundingBox() { return kdtree ? kdtree->box : AABB(); };

private:
	bool initFromScene(const aiScene* pScene, const std::string& Filename);
	void initMesh(unsigned int index, const aiMesh* paiMesh);
	void initMaterial(unsigned int index, const aiMaterial *material);

	void initTriangles();
	void clear();

	void initKDTree();

	std::vector<Mesh> meshes;
	std::vector<Material> materials;

	std::vector<Triangle*> tris;

	KDNode *kdtree;
};