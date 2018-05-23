#ifndef MESH_TOPOLOGY_H
#define MESH_TOPOLOGY_H
#include <map>
#include <algorithm>
#include <vector>
#include "vector.h"
#include "optimizedModel.h"

struct TopologyData
{
	float *pTopoVertices;
	int vertexSize;

	unsigned int *pTopoIndices;
	int indexSize;

	void AllocVertex(int size)
	{
		pTopoVertices = new float[size];
		vertexSize = size;
	}

	void AllocIndex(int size)
	{
		pTopoIndices = new unsigned int[size];
		indexSize = size;
	}

	void ReleaseVertex()
	{
		if (pTopoVertices)
		{
			delete[]pTopoVertices;
			pTopoVertices = NULL;
			vertexSize = 0;
		}
	}

	void ReleaseIndex()
	{
		if (pTopoIndices)
		{
			delete[]pTopoIndices;
			pTopoIndices = NULL;
			indexSize = 0;
		}
	}

	void SetOptimizedVolume(OptimizedVolume *ov)
	{
		AllocIndex(ov->faces.size() * 3);
		AllocVertex(ov->points.size() * 3);

		for (unsigned int i = 0; i < ov->faces.size(); i++)
		{
			int idx = i * 3;
			pTopoIndices[idx] = ov->faces[i].index[0];
			pTopoIndices[idx + 1] = ov->faces[i].index[1];
			pTopoIndices[idx + 2] = ov->faces[i].index[2];
		}

		for (unsigned int i = 0; i < ov->points.size(); i++)
		{
			int idx = i * 3;
			pTopoVertices[idx] = ov->points[i].p.x;
			pTopoVertices[idx + 1] = ov->points[i].p.y;
			pTopoVertices[idx + 2] = ov->points[i].p.z;
		}
	}
};

class FaceGroup
{
private:
	std::vector<unsigned int> mSharedFaces;

public:
	FaceGroup()
	{
	}

	~FaceGroup()
	{
	}

	void AddFace(unsigned int faceIndex)
	{
		mSharedFaces.push_back(faceIndex);
	}

	Vector3 CalcNormal(float *refFaceNormals)
	{
		Vector3 normal(0, 0, 0);
		for (int i = 0; i < mSharedFaces.size(); i++)
		{
			unsigned int idx = (mSharedFaces[i] * 3);
			normal+= Vector3(refFaceNormals[idx], refFaceNormals[idx + 1], refFaceNormals[idx + 2]);
		}
		normal.NormalizeSafe();
		return normal;
	}

	std::vector<unsigned int> *GetFaceListPtr()
	{
		return &mSharedFaces;
	}
};

struct QuadFace
{
	unsigned int A, B, C, D;
	int FaceA, FaceB;
	QuadFace(unsigned int a, unsigned int b, unsigned int c, unsigned int d, int faceA, int faceB)
	{
		A = a;
		B = b;
		C = c;
		D = d;
		FaceA = faceA;
		FaceB = faceB;
	}
};

struct FaceData
{
	int Index;
	unsigned int A;
	unsigned int B;
	unsigned int C;

	FaceData(unsigned int a=0, unsigned int b=0, unsigned int c=0, int face=-1)
	{
		A = a;
		B = b;
		C = c;
		Index = face;
	}
};

TopologyData BuildVertexTopology(float *pVertices, int vertexSize);
std::vector<FaceGroup> *MakeNeighborFaceOfVertex(TopologyData &topoData);
int FindNeighborFaceOfEdge(const std::vector<FaceGroup> &neighborFaceList, int vertexIndex0, int vertexIndex1);
void SeparateTriAndQuadFace(TopologyData *pTopoData, std::vector<QuadFace> *pQuadList, std::vector<FaceData> *pTriList);

#endif//MESH_TOPOLOGY_H