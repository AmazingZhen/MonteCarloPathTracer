#include "Model.h"


bool Model::loadModel(const std::string & filename)
{
	// Release the previously loaded Model (if it exists)
	clear();

	bool ret = false;
	Assimp::Importer importer;

	const aiScene* pScene = importer.ReadFile(filename.c_str(), aiProcess_Triangulate | aiProcess_GenSmoothNormals);

	if (pScene) {
		ret = initFromScene(pScene, filename);
	}
	else {
		printf("Error parsing '%s': '%s'\n", filename.c_str(), importer.GetErrorString());
	}

	return ret;
}

bool Model::intersect(const Ray & r, IntersectionInfo & info)
{
	float t, u, v, tmin = FLT_MAX;
	return KDNode::intersect(kdtree, r, tmin, info);

	/*
	for (Triangle* tri : tris) if (tri->intersect(r, t, u, v) && t < tmin) {
		tmin = t;

		tri->getSurfaceProperties(r, t, u, v, info);
	}

	return info.hit;
	*/
}

bool Model::initFromScene(const aiScene * pScene, const std::string & Filename)
{
	meshes.resize(pScene->mNumMeshes);
	materials.resize(pScene->mNumMaterials);

	for (unsigned int i = 0; i < materials.size(); i++) {
		const aiMaterial* paiMaterial = pScene->mMaterials[i];
		initMaterial(i, paiMaterial);
	}

	// Initialize the Modeles in the scene one by one
	for (unsigned int i = 0; i < meshes.size(); i++) {
		const aiMesh* paiMesh = pScene->mMeshes[i];
		initMesh(i, paiMesh);
	}

	initTriangles();

	initKDTree();

	return true;
}

void Model::initMesh(unsigned int index, const aiMesh *paiMesh)
{
	//if (index != 3) {
	//	return;
	//}

	meshes[index].materialIndex = paiMesh->mMaterialIndex;

	std::vector<glm::vec3> &vertices_ = meshes[index].vertices;
	std::vector<unsigned int> &indices_ = meshes[index].indices;
	std::vector<glm::vec3> &normals_ = meshes[index].normals;

	const aiVector3D Zero3D(0.0f, 0.0f, 0.0f);

	for (unsigned int i = 0; i < paiMesh->mNumVertices; i++) {
		const aiVector3D* pPos = &(paiMesh->mVertices[i]);
		const aiVector3D* pNormal = &(paiMesh->mNormals[i]);
		// const aiVector3D* pTexCoord = paiMesh->HasTextureCoords(0) ? &(paiMesh->mTextureCoords[0][i]) : &Zero3D;

		vertices_.push_back(glm::vec3(pPos->x, pPos->y, pPos->z));
		normals_.push_back(glm::vec3(pNormal->x, pNormal->y, pNormal->z));
	}

	for (unsigned int i = 0; i < paiMesh->mNumFaces; i++) {
		const aiFace& face = paiMesh->mFaces[i];
		assert(face.mNumIndices == 3);
		indices_.push_back(face.mIndices[0]);
		indices_.push_back(face.mIndices[1]);
		indices_.push_back(face.mIndices[2]);
	}
}

void Model::initMaterial(unsigned int index, const aiMaterial *paiMaterial)
{
	Material *mater = &(materials[index]);

	aiString mname;
	paiMaterial->Get(AI_MATKEY_NAME, mname);
	if (mname.length > 0)
		mater->name = mname.C_Str();

	int shadingModel;
	paiMaterial->Get(AI_MATKEY_SHADING_MODEL, shadingModel);

	if (shadingModel != aiShadingMode_Phong && shadingModel != aiShadingMode_Gouraud)
	{
		printf("This Model's shading model is not implemented in this loader, setting to default paiMaterial\n");
		mater->name = "DefaultMaterial";
	}
	else
	{
		aiColor3D dif(0.f, 0.f, 0.f);
		aiColor3D amb(0.f, 0.f, 0.f);
		aiColor3D spec(0.f, 0.f, 0.f);
		aiColor3D emiss(0.f, 0.f, 0.f);

		float shine = 0.f;
		float refractIndex = 0.f;

		aiColor3D trans(0.f, 0.f, 0.f);

		paiMaterial->Get(AI_MATKEY_COLOR_AMBIENT, amb);
		paiMaterial->Get(AI_MATKEY_COLOR_DIFFUSE, dif); //->Get(<paiMaterial-key>,<where-to-store>))
		paiMaterial->Get(AI_MATKEY_COLOR_SPECULAR, spec);
		paiMaterial->Get(AI_MATKEY_COLOR_EMISSIVE, emiss);

		paiMaterial->Get(AI_MATKEY_SHININESS, shine);
		paiMaterial->Get(AI_MATKEY_REFRACTI, refractIndex);

		paiMaterial->Get(AI_MATKEY_COLOR_TRANSPARENT, trans);

		mater->ka = glm::vec3(amb.r, amb.g, amb.b);
		mater->kd = glm::vec3(dif.r, dif.g, dif.b);
		mater->ks = glm::vec3(spec.r, spec.g, spec.b);
		mater->ke = glm::vec3(emiss.r, emiss.g, emiss.b);

		mater->ns = shine;
		mater->ni = refractIndex;

		mater->tf = glm::vec3(trans.r, trans.g, trans.b);
		/*
		mater->ambient *= .2f;
		if (mater->shininess == 0.0) {
			mater->shininess = 30;
		}
		*/
	}

	mater->print();
}

void Model::initTriangles()
{
	for (Mesh &m : meshes) {
		for (int i = 0; i < m.indices.size(); i += 3) {
			unsigned int index[3];

			index[0] = m.indices[i];
			index[1] = m.indices[i + 1];
			index[2] = m.indices[i + 2];

			Triangle *tri = new Triangle(index, &m, &(materials[m.materialIndex]));

			tris.push_back(tri);
		}
	}
}

void Model::clear()
{
	meshes.clear();
	materials.clear();
}

void Model::initKDTree()
{
	kdtree = KDNode::build(tris, 0);
}

