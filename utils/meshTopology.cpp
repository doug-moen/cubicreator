#include "meshTopology.h"
#include "util.h"

TopologyData BuildVertexTopology(float *pVertices, int vertexSize)
{
	std::vector<float> vertexList;
	std::vector<unsigned int> indexList;
	std::map<int, unsigned int> map;
	//std::map<std::string, unsigned int> map;
	int MELD_DIST = 30;

	for (unsigned int i = 0; i < vertexSize / 3; i++)
	{
		unsigned int idx = i * 3;
		float a = pVertices[idx];
		float b = pVertices[idx + 1];
		float c = pVertices[idx + 2];
		Vector3 p(a, b, c);

		int x = (int)(p.X * 1000);
		int y = (int)(p.Y * 1000);
		int z = (int)(p.Z * 1000);
		int key = ((x + MELD_DIST / 2) / MELD_DIST) ^ (((y + MELD_DIST / 2) / MELD_DIST) << 10) ^ (((z + MELD_DIST / 2) / MELD_DIST) << 20);
		unsigned int index;
		std::map<int, unsigned int>::iterator it = map.find(key);
		//std::map<std::string, unsigned int>::iterator it = map.find(key);
		if (it != map.end())
		{
			//map.TryGetValue(key, out index);
			index = it->second;
		}
		else
		{
			//add here
			index = ((unsigned int)vertexList.size() / 3);
			vertexList.push_back(a);
			vertexList.push_back(b);
			vertexList.push_back(c);
			map.insert({ key, index });
		}
		indexList.push_back(index);
	}

	map.clear();
	TopologyData topoData;
	topoData.AllocVertex(vertexList.size());
	std::copy(vertexList.begin(), vertexList.end(), topoData.pTopoVertices);

	topoData.AllocIndex(indexList.size());
	std::copy(indexList.begin(), indexList.end(), topoData.pTopoIndices);
	return topoData;
}

std::vector<FaceGroup> *MakeNeighborFaceOfVertex(TopologyData &topoData)
{
	int totalVertex = topoData.vertexSize / 3;

	std::vector<FaceGroup> *pNeighborList = new std::vector<FaceGroup>;
	std::vector<FaceGroup> &neighborList = *pNeighborList;
	//pNeighborList->resize(totalVertex);
	for (int i = 0; i < totalVertex; i++)
		pNeighborList->push_back(FaceGroup());

	for (unsigned int face = 0; face < topoData.indexSize / 3; face++)
	{
		unsigned int idx = face * 3;
		int a = (int)topoData.pTopoIndices[idx	  ];
		int b = (int)topoData.pTopoIndices[idx + 1];
		int c = (int)topoData.pTopoIndices[idx + 2];

		neighborList[a].AddFace(face);
		neighborList[b].AddFace(face);
		neighborList[c].AddFace(face);
	}
	return pNeighborList;
}

int FindNeighborFaceOfEdge(std::vector<FaceGroup> &neighborFaceList, int vertexIndex0, int vertexIndex1)
{
	std::vector<unsigned int> &faceList = *neighborFaceList[vertexIndex0].GetFaceListPtr();

	for (int i = 0; i< faceList.size(); i++)
	{
		unsigned int faceIndex = faceList[i];

		std::vector<unsigned int> &faceList1 = *neighborFaceList[vertexIndex1].GetFaceListPtr();
		for (int k = 0; k<faceList1.size(); k++)
		{
			if (faceIndex == faceList1[k])
				return (int)faceIndex;
		}
	}
	return -1;
}

FaceData GetFaceData(unsigned int *pIndices, unsigned int face)
{
	FaceData faceData;
	faceData.Index = face;
	unsigned int idx = face * 3;
	faceData.A = pIndices[idx];
	faceData.B = pIndices[idx + 1];
	faceData.C = pIndices[idx + 2];
	return faceData;
}

float FindOppositeAngle(float *topoVertices, unsigned int a, unsigned int b, FaceData faceData, unsigned int *pOutD, Vector3 srcNormal)
{
	Vector3 p0 = Vector3(topoVertices[faceData.A * 3], topoVertices[faceData.A * 3 + 1], topoVertices[faceData.A * 3 + 2]);
	Vector3 p1 = Vector3(topoVertices[faceData.B * 3], topoVertices[faceData.B * 3 + 1], topoVertices[faceData.B * 3 + 2]);
	Vector3 p2 = Vector3(topoVertices[faceData.C * 3], topoVertices[faceData.C * 3 + 1], topoVertices[faceData.C * 3 + 2]);
	Vector3 normal = calcNormal(p0, p1, p2);
	*pOutD = 0;
	float deg = 0;

	//if (srcNormal != normal)
	if (abs(calcAngle(srcNormal, normal))>0.2f)
		return 0;

	if (a != faceData.A && b != faceData.A)
	{
		deg = calcSignedAngle(p1 - p0, p2 -p0, normal);
		*pOutD = faceData.A;
	}

	if (a != faceData.B && b != faceData.B)
	{
		deg = calcSignedAngle(p0 -p1, p2-p1, normal);
		*pOutD = faceData.B;
	}

	if (a != faceData.C && b != faceData.C)
	{
		deg = calcSignedAngle(p0 - p2, p1 - p2, normal);
		*pOutD = faceData.C;
	}
	return deg;
}

void SeparateTriAndQuadFace(TopologyData *pTopoData, std::vector<QuadFace> *pQuadList, std::vector<FaceData> *pTriList)
{
	std::vector<QuadFace> &quadList = *pQuadList;
	std::vector<FaceData> &triList = *pTriList;
		
	TopologyData &topoData = *pTopoData;// BuildVertexTopology(pVertices, vertexSize);//cstyle 이미 만들어져 있다.
	std::vector<FaceGroup> *pNeighborFaceList = MakeNeighborFaceOfVertex(topoData);
	std::vector<FaceGroup> &neighborFaceList = *pNeighborFaceList;

	//1. find neighbor faces
	//2. distinguish face can be a quad with neighbors.
	//3. if a face is possible to be a quad, then put it to a quad list and delete the faces.
	//4. when applied step 3 to all faces, then rest faces put to the TriFace list.
	std::vector<bool> checkList;
	for (unsigned int i = 0; i < topoData.indexSize / 3; i++)
		checkList.push_back(true);
	//make quad list
	for (int i = 0; i < checkList.size(); i++)
	{
		int curFaceIndex = i;
		//check if the face is already processed.
		if (!checkList[curFaceIndex])
			continue;

		unsigned int idx = (unsigned int)i * 3;
		unsigned int viA = topoData.pTopoIndices[idx	];
		unsigned int viB = topoData.pTopoIndices[idx + 1];
		unsigned int viC = topoData.pTopoIndices[idx + 2];

		unsigned int a = viA * 3;
		unsigned int b = viB * 3;
		unsigned int c = viC * 3;

		//Point3로 변경
		Vector3 va = Vector3(topoData.pTopoVertices[a], topoData.pTopoVertices[a + 1], topoData.pTopoVertices[a + 2]);
		Vector3 vb = Vector3(topoData.pTopoVertices[b], topoData.pTopoVertices[b + 1], topoData.pTopoVertices[b + 2]);
		Vector3 vc = Vector3(topoData.pTopoVertices[c], topoData.pTopoVertices[c + 1], topoData.pTopoVertices[c + 2]);
		Vector3 normal = calcNormal(va, vb, vc);//double로 계산?
		
		//p0 - p1, p2에 걸친면 중복되지 않는 점을 찾고 거기에 각도를 구해서 80도 이상이면 사각형을 만들수 있는 면이 된다.
		//va에 상대되는  p를 찾는다.
		//px의 각도를 구해서 80도가 넘으면 사각형만들수 있다.
		//사각형이 가능하면 curFaceIndex 와 faceIndexBC를 checkList에서 지운다.
		float minAngle = 80;
		float maxAngle = 100;
		//a 
		float angle = calcSignedAngle(vb - va, vc - va, normal);
		if (abs(angle)>minAngle)
		{
			int faceIndexBC = FindNeighborFaceOfEdge(neighborFaceList, (int)viB, (int)viC);
			if (faceIndexBC != curFaceIndex && checkList[faceIndexBC])
			{
				FaceData faceData = GetFaceData(topoData.pTopoIndices, (unsigned int)faceIndexBC);
				unsigned int d;
				if (abs(FindOppositeAngle(topoData.pTopoVertices, viB, viC, faceData, &d, normal)) > minAngle)
				{					
					quadList.push_back(QuadFace(viA, viB, d, viC, curFaceIndex, faceIndexBC));	
					checkList[curFaceIndex] = false;
					checkList[faceIndexBC] = false;
				}
			}
		}
		//b
		angle = calcSignedAngle(va -vb, vc-vb, normal);
		if (abs(angle) > minAngle)
		{
			int faceIndexAC = FindNeighborFaceOfEdge(neighborFaceList, (int)viA, (int)viC);
			if (faceIndexAC != curFaceIndex && checkList[faceIndexAC])
			{
				FaceData faceData = GetFaceData(topoData.pTopoIndices, (unsigned int)faceIndexAC);
				unsigned int d;
				if (abs(FindOppositeAngle(topoData.pTopoVertices, viA, viC, faceData, &d, normal)) > minAngle)
				{					
					quadList.push_back(QuadFace(viB, viC, d, viA, curFaceIndex, faceIndexAC));
					checkList[curFaceIndex] = false;
					checkList[faceIndexAC] = false;
				}
			}
		}
		//c
		angle = calcSignedAngle(va - vc, vb-vc, normal);
		if (abs(angle) > minAngle)
		{
			int faceIndexAB = FindNeighborFaceOfEdge(neighborFaceList, (int)viA, (int)viB);
			if (faceIndexAB != curFaceIndex && checkList[faceIndexAB])
			{
				FaceData faceData = GetFaceData(topoData.pTopoIndices, (unsigned int)faceIndexAB);
				unsigned int d;
				if (abs(FindOppositeAngle(topoData.pTopoVertices, viA, viB, faceData, &d, normal)) > minAngle)
				{					
					quadList.push_back(QuadFace(viC, viA, d, viB, curFaceIndex, faceIndexAB));
					checkList[curFaceIndex] = false;
					checkList[faceIndexAB] = false;
				}
			}
		}
	}

	//make triangle list
	for (int i = 0; i < checkList.size(); i++)
	{
		if (checkList[i])		
			triList.push_back(FaceData(topoData.pTopoIndices[i * 3], topoData.pTopoIndices[i * 3 + 1], topoData.pTopoIndices[i * 3 + 2], i));					
	}
}