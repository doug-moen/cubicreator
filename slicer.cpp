/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <stdio.h>

#include "utils/gettime.h"
#include "utils/logoutput.h"

#include "sliceDataStorage.h"
#include "slicer.h"
#include "polygonOptimizer.h"
#include "settings.h"

#include "layerPart.h"
#include "support.h"
#include "utils/meshTopology.h"

void SlicerLayer::makePolygons(OptimizedVolume* ov, bool keepNoneClosed, bool extensiveStitching)
{
    for(unsigned int startSegment=0; startSegment < segmentList.size(); startSegment++)
    {
        if (segmentList[startSegment].addedToPolygon)
            continue;
        
        Polygon poly;
        poly.add(segmentList[startSegment].start);
        
        unsigned int segmentIndex = startSegment;
        bool canClose;
        while(true)
        {
            canClose = false;
            segmentList[segmentIndex].addedToPolygon = true;
            Point p0 = segmentList[segmentIndex].end;
            poly.add(p0);
            int nextIndex = -1;
            OptimizedFace* face = &ov->faces[segmentList[segmentIndex].faceIndex];
            for(unsigned int i=0;i<3;i++)
            {
                if (face->touching[i] > -1 && faceToSegmentIndex.find(face->touching[i]) != faceToSegmentIndex.end())
                {
     //               Point p1 = segmentList[faceToSegmentIndex[face->touching[i]]].start;
     //               Point diff = p0 - p1;
     //               if (shorterThen(diff, 10))//org
					////if (shorterThen(diff, MINIMAL_POLYGON_LENGTH))
     //               {
     //                   if (faceToSegmentIndex[face->touching[i]] == (int)startSegment)
     //                       canClose = true;
     //                   if (segmentList[faceToSegmentIndex[face->touching[i]]].addedToPolygon)
     //                       continue;
     //                   nextIndex = faceToSegmentIndex[face->touching[i]];
     //               }

					int segIndex = faceToSegmentIndex[face->touching[i]];
					if (segmentList[segIndex].addedToPolygon || (segmentList[segmentIndex].faceIndex1 != -1 && segmentList[segmentIndex].faceIndex1 == face->touching[i]))
						continue;
					Point p1 = segmentList[segIndex].start;
					Point diff = p0 - p1;
					if (shorterThen(diff, 10))//org
						//if (shorterThen(diff, MINIMAL_POLYGON_LENGTH))
					{
						//if (faceToSegmentIndex[face->touching[i]] == (int)startSegment)
						if (segIndex == (int)startSegment)
							canClose = true;
						//if (segmentList[faceToSegmentIndex[face->touching[i]]].addedToPolygon)
						//    continue;
						nextIndex = segIndex;//faceToSegmentIndex[face->touching[i]];
					}
                }
            }
            if (nextIndex == -1)
                break;
            segmentIndex = nextIndex;
        }
        if (canClose)
            polygonList.add(poly);
        else
            openPolygonList.add(poly);
    }
    //Clear the segmentList to save memory, it is no longer needed after this point.
    segmentList.clear();

    //Connecting polygons that are not closed yet, as models are not always perfect manifold we need to join some stuff up to get proper polygons
    //First link up polygon ends that are within 2 microns.
    for(unsigned int i=0;i<openPolygonList.size();i++)
    {
        if (openPolygonList[i].size() < 1) continue;
        for(unsigned int j=0;j<openPolygonList.size();j++)
        {
            if (openPolygonList[j].size() < 1) continue;
            
            Point diff = openPolygonList[i][openPolygonList[i].size()-1] - openPolygonList[j][0];
            int64_t distSquared = vSize2(diff);

            if (distSquared < 2 * 2)
            {
                if (i == j)
                {
                    polygonList.add(openPolygonList[i]);
                    openPolygonList[i].clear();
                    break;
                }else{
                    for(unsigned int n=0; n<openPolygonList[j].size(); n++)
                        openPolygonList[i].add(openPolygonList[j][n]);

                    openPolygonList[j].clear();
                }
            }
        }
    }

    //Next link up all the missing ends, closing up the smallest gaps first. This is an inefficient implementation which can run in O(n*n*n) time.
    while(1)
    {
        int64_t bestScore = 10000 * 10000;
        unsigned int bestA = -1;
        unsigned int bestB = -1;
        bool reversed = false;
        for(unsigned int i=0;i<openPolygonList.size();i++)
        {
            if (openPolygonList[i].size() < 1) continue;
            for(unsigned int j=0;j<openPolygonList.size();j++)
            {
                if (openPolygonList[j].size() < 1) continue;
                
                Point diff = openPolygonList[i][openPolygonList[i].size()-1] - openPolygonList[j][0];
                int64_t distSquared = vSize2(diff);
                if (distSquared < bestScore)
                {
                    bestScore = distSquared;
                    bestA = i;
                    bestB = j;
                    reversed = false;
                }

                if (i != j)
                {
                    Point diff = openPolygonList[i][openPolygonList[i].size()-1] - openPolygonList[j][openPolygonList[j].size()-1];
                    int64_t distSquared = vSize2(diff);
                    if (distSquared < bestScore)
                    {
                        bestScore = distSquared;
                        bestA = i;
                        bestB = j;
                        reversed = true;
                    }
                }
            }
        }
        
        if (bestScore >= 10000 * 10000)
            break;
        
        if (bestA == bestB)
        {
            polygonList.add(openPolygonList[bestA]);
            openPolygonList[bestA].clear();
        }else{
            if (reversed)
            {
                if (openPolygonList[bestA].polygonLength() > openPolygonList[bestB].polygonLength())
                {
                    for(unsigned int n=openPolygonList[bestB].size()-1; int(n)>=0; n--)
                        openPolygonList[bestA].add(openPolygonList[bestB][n]);
                    openPolygonList[bestB].clear();
                }else{
                    for(unsigned int n=openPolygonList[bestA].size()-1; int(n)>=0; n--)
                        openPolygonList[bestB].add(openPolygonList[bestA][n]);
                    openPolygonList[bestA].clear();
                }
            }else{
                for(unsigned int n=0; n<openPolygonList[bestB].size(); n++)
                    openPolygonList[bestA].add(openPolygonList[bestB][n]);
                openPolygonList[bestB].clear();
            }
        }
    }

    if (extensiveStitching)
    {
        //For extensive stitching find 2 open polygons that are touching 2 closed polygons.
        // Then find the sortest path over this polygon that can be used to connect the open polygons,
        // And generate a path over this shortest bit to link up the 2 open polygons.
        // (If these 2 open polygons are the same polygon, then the final result is a closed polyon)
        
        while(1)
        {
            unsigned int bestA = -1;
            unsigned int bestB = -1;
            gapCloserResult bestResult;
            bestResult.len = LLONG_MAX;
            bestResult.polygonIdx = -1;
            bestResult.pointIdxA = -1;
            bestResult.pointIdxB = -1;
            
            for(unsigned int i=0; i<openPolygonList.size(); i++)
            {
                if (openPolygonList[i].size() < 1) continue;
                
                {
                    gapCloserResult res = findPolygonGapCloser(openPolygonList[i][0], openPolygonList[i][openPolygonList[i].size()-1]);
                    if (res.len > 0 && res.len < bestResult.len)
                    {
                        bestA = i;
                        bestB = i;
                        bestResult = res;
                    }
                }

                for(unsigned int j=0; j<openPolygonList.size(); j++)
                {
                    if (openPolygonList[j].size() < 1 || i == j) continue;
                    
                    gapCloserResult res = findPolygonGapCloser(openPolygonList[i][0], openPolygonList[j][openPolygonList[j].size()-1]);
                    if (res.len > 0 && res.len < bestResult.len)
                    {
                        bestA = i;
                        bestB = j;
                        bestResult = res;
                    }
                }
            }
            
            if (bestResult.len < LLONG_MAX)
            {
                if (bestA == bestB)
                {
                    if (bestResult.pointIdxA == bestResult.pointIdxB)
                    {
                        polygonList.add(openPolygonList[bestA]);
                        openPolygonList[bestA].clear();
                    }
                    else if (bestResult.AtoB)
                    {
                        PolygonRef poly = polygonList.newPoly();
                        for(unsigned int j = bestResult.pointIdxA; j != bestResult.pointIdxB; j = (j + 1) % polygonList[bestResult.polygonIdx].size())
                            poly.add(polygonList[bestResult.polygonIdx][j]);
                        for(unsigned int j = openPolygonList[bestA].size() - 1; int(j) >= 0; j--)
                            poly.add(openPolygonList[bestA][j]);
                        openPolygonList[bestA].clear();
                    }
                    else
                    {
                        unsigned int n = polygonList.size();
                        polygonList.add(openPolygonList[bestA]);
                        for(unsigned int j = bestResult.pointIdxB; j != bestResult.pointIdxA; j = (j + 1) % polygonList[bestResult.polygonIdx].size())
                            polygonList[n].add(polygonList[bestResult.polygonIdx][j]);
                        openPolygonList[bestA].clear();
                    }
                }
                else
                {
                    if (bestResult.pointIdxA == bestResult.pointIdxB)
                    {
                        for(unsigned int n=0; n<openPolygonList[bestA].size(); n++)
                            openPolygonList[bestB].add(openPolygonList[bestA][n]);
                        openPolygonList[bestA].clear();
                    }
                    else if (bestResult.AtoB)
                    {
                        Polygon poly;
                        for(unsigned int n = bestResult.pointIdxA; n != bestResult.pointIdxB; n = (n + 1) % polygonList[bestResult.polygonIdx].size())
                            poly.add(polygonList[bestResult.polygonIdx][n]);
                        for(unsigned int n=poly.size()-1;int(n) >= 0; n--)
                            openPolygonList[bestB].add(poly[n]);
                        for(unsigned int n=0; n<openPolygonList[bestA].size(); n++)
                            openPolygonList[bestB].add(openPolygonList[bestA][n]);
                        openPolygonList[bestA].clear();
                    }
                    else
                    {
                        for(unsigned int n = bestResult.pointIdxB; n != bestResult.pointIdxA; n = (n + 1) % polygonList[bestResult.polygonIdx].size())
                            openPolygonList[bestB].add(polygonList[bestResult.polygonIdx][n]);
                        for(unsigned int n = openPolygonList[bestA].size() - 1; int(n) >= 0; n--)
                            openPolygonList[bestB].add(openPolygonList[bestA][n]);
                        openPolygonList[bestA].clear();
                    }
                }
            }
            else
            {
                break;
            }
        }
    }

    /*
    int q=0;
    for(unsigned int i=0;i<openPolygonList.size();i++)
    {
        if (openPolygonList[i].size() < 2) continue;
        if (!q) log("***\n");
        log("S: %f %f\n", float(openPolygonList[i][0].X), float(openPolygonList[i][0].Y));
        log("E: %f %f\n", float(openPolygonList[i][openPolygonList[i].size()-1].X), float(openPolygonList[i][openPolygonList[i].size()-1].Y));
        q = 1;
    }
    */
    //if (q) exit(1);

    if (keepNoneClosed)
    {
        for(unsigned int n=0; n<openPolygonList.size(); n++)
        {
            if (openPolygonList[n].size() > 0)
                polygonList.add(openPolygonList[n]);
        }
    }
    //Clear the openPolygonList to save memory, the only reason to keep it after this is for debugging.
    //openPolygonList.clear();

    //Remove all the tiny polygons, or polygons that are not closed. As they do not contribute to the actual print.
    int snapDistance = 1000;
    for(unsigned int i=0;i<polygonList.size();i++)
    {
        int length = 0;
        
        for(unsigned int n=1; n<polygonList[i].size(); n++)
        {
            length += vSize(polygonList[i][n] - polygonList[i][n-1]);
            if (length > snapDistance)
                break;
        }
        if (length < snapDistance)
        {
            polygonList.remove(i);
            i--;
        }
    }

    //Finally optimize all the polygons. Every point removed saves time in the long run.
    optimizePolygons(polygonList);
	//polygonList.simplify();
	//polygonList.removeDegenerateVerts(); // remove verts connected to overlapping line segments
}


Slicer::Slicer(OptimizedVolume* ov, int32_t initial, int32_t thickness, bool keepNoneClosed, bool extensiveStitching)
{
    modelSize = ov->model->modelSize;
    modelMin = ov->model->vMin;
    
    int layerCount = (modelSize.z - initial) / thickness + 1;
    //log("Layer count: %i\n", layerCount);
    layers.resize(layerCount);
    
    for(int32_t layerNr = 0; layerNr < layerCount; layerNr++)
    {
        layers[layerNr].z = initial + thickness * layerNr;
    }
    
    for(unsigned int i=0; i<ov->faces.size(); i++)
    {
        Point3 p0 = ov->points[ov->faces[i].index[0]].p;
        Point3 p1 = ov->points[ov->faces[i].index[1]].p;
        Point3 p2 = ov->points[ov->faces[i].index[2]].p;
        int32_t minZ = p0.z;
        int32_t maxZ = p0.z;
        if (p1.z < minZ) minZ = p1.z;
        if (p2.z < minZ) minZ = p2.z;
        if (p1.z > maxZ) maxZ = p1.z;
        if (p2.z > maxZ) maxZ = p2.z;
        
        for(int32_t layerNr = (minZ - initial) / thickness; layerNr <= (maxZ - initial) / thickness; layerNr++)
        {
            int32_t z = layerNr * thickness + initial;
            if (z < minZ) continue;
            if (layerNr < 0) continue;
            
            SlicerSegment s;
            if (p0.z < z && p1.z >= z && p2.z >= z)
                s = project2D(p0, p2, p1, z);
            else if (p0.z > z && p1.z < z && p2.z < z)
                s = project2D(p0, p1, p2, z);

            else if (p1.z < z && p0.z >= z && p2.z >= z)
                s = project2D(p1, p0, p2, z);
            else if (p1.z > z && p0.z < z && p2.z < z)
                s = project2D(p1, p2, p0, z);

            else if (p2.z < z && p1.z >= z && p0.z >= z)
                s = project2D(p2, p1, p0, z);
            else if (p2.z > z && p1.z < z && p0.z < z)
                s = project2D(p2, p0, p1, z);
            else
            {
                //Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
                //  on the slice would create two segments
                continue;
            }
            layers[layerNr].faceToSegmentIndex[i] = layers[layerNr].segmentList.size();
            s.faceIndex = i;
            s.addedToPolygon = false;
            layers[layerNr].segmentList.push_back(s);
        }
    }
    
    for(unsigned int layerNr=0; layerNr<layers.size(); layerNr++)
    {
        layers[layerNr].makePolygons(ov, keepNoneClosed, extensiveStitching);
    }
}

//Slicer::Slicer(OptimizedVolume* ov, int32_t initial, int32_t thickness, bool keepNoneClosed, bool extensiveStitching)
//{
//	modelSize = ov->model->modelSize;
//	modelMin = ov->model->vMin;
//
//	int layerCount = (modelSize.z - initial) / thickness + 1;
//	log("Layer count: %i\n", layerCount);
//	layers.resize(layerCount);
//
//	for (int32_t layerNr = 0; layerNr < layerCount; layerNr++)
//		layers[layerNr].z = initial + thickness * layerNr;
//
//	//cstyle make topoData from OptimizedVolume
//	TopologyData topoData;
//	topoData.SetOptimizedVolume(ov);
//	std::vector<QuadFace> quadList;
//	std::vector<FaceData> triList;
//	SeparateTriAndQuadFace(&topoData, &quadList, &triList);
//
//	//slice triangle face
//	for (unsigned int i = 0; i<triList.size(); i++)
//	{
//		int idxA = triList[i].A * 3;
//		int idxB = triList[i].B * 3;
//		int idxC = triList[i].C * 3;
//		Point3 p0 = Point3(topoData.pTopoVertices[idxA], topoData.pTopoVertices[idxA + 1], topoData.pTopoVertices[idxA + 2]);
//		Point3 p1 = Point3(topoData.pTopoVertices[idxB], topoData.pTopoVertices[idxB + 1], topoData.pTopoVertices[idxB + 2]);
//		Point3 p2 = Point3(topoData.pTopoVertices[idxC], topoData.pTopoVertices[idxC + 1], topoData.pTopoVertices[idxC + 2]);
//
//		int32_t minZ = p0.z;
//		int32_t maxZ = p0.z;
//
//		if (p1.z < minZ)
//			minZ = p1.z;
//		if (p2.z < minZ)
//			minZ = p2.z;
//		if (p1.z > maxZ)
//			maxZ = p1.z;
//		if (p2.z > maxZ)
//			maxZ = p2.z;
//
//		for (int32_t layerNr = (minZ - initial) / thickness; layerNr <= (maxZ - initial) / thickness; layerNr++)
//		{
//			int32_t z = layerNr * thickness + initial;
//			if (z < minZ)
//				continue;
//			if (layerNr < 0)
//				continue;
//
//			SlicerSegment s;
//			if (p0.z < z && p1.z >= z && p2.z >= z)
//				s = project2D(p0, p2, p1, z);
//			else if (p0.z > z && p1.z < z && p2.z < z)
//				s = project2D(p0, p1, p2, z);
//
//			else if (p1.z < z && p0.z >= z && p2.z >= z)
//				s = project2D(p1, p0, p2, z);
//			else if (p1.z > z && p0.z < z && p2.z < z)
//				s = project2D(p1, p2, p0, z);
//
//			else if (p2.z < z && p1.z >= z && p0.z >= z)
//				s = project2D(p2, p1, p0, z);
//			else if (p2.z > z && p1.z < z && p0.z < z)
//				s = project2D(p2, p0, p1, z);
//			else
//			{
//				//Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
//				//  on the slice would create two segments
//				continue;
//			}
//
//			int face = triList[i].Index;
//			layers[layerNr].faceToSegmentIndex[face] = layers[layerNr].segmentList.size();
//			s.faceIndex = face;
//			s.addedToPolygon = false;
//			layers[layerNr].segmentList.push_back(s);
//		}
//	}
//
//	//slice quad face
//	for (unsigned int i = 0; i<quadList.size(); i++)
//	{
//		int idxA = quadList[i].A * 3;
//		int idxB = quadList[i].B * 3;
//		int idxC = quadList[i].C * 3;
//		int idxD = quadList[i].D * 3;
//		int faceA = quadList[i].FaceA;
//		int faceB = quadList[i].FaceB;
//
//		Point3 p0 = Point3(topoData.pTopoVertices[idxA], topoData.pTopoVertices[idxA + 1], topoData.pTopoVertices[idxA + 2]);
//		Point3 p1 = Point3(topoData.pTopoVertices[idxB], topoData.pTopoVertices[idxB + 1], topoData.pTopoVertices[idxB + 2]);
//		Point3 p2 = Point3(topoData.pTopoVertices[idxC], topoData.pTopoVertices[idxC + 1], topoData.pTopoVertices[idxC + 2]);
//		Point3 p3 = Point3(topoData.pTopoVertices[idxD], topoData.pTopoVertices[idxD + 1], topoData.pTopoVertices[idxD + 2]);
//
//		int32_t minZ = p0.z;
//		int32_t maxZ = p0.z;
//
//		if (p1.z < minZ)
//			minZ = p1.z;
//		if (p2.z < minZ)
//			minZ = p2.z;
//		if (p3.z < minZ)
//			minZ = p3.z;
//
//		if (p1.z > maxZ)
//			maxZ = p1.z;
//		if (p2.z > maxZ)
//			maxZ = p2.z;
//		if (p3.z > maxZ)
//			maxZ = p3.z;
//
//		for (int32_t layerNr = (minZ - initial) / thickness; layerNr <= (maxZ - initial) / thickness; layerNr++)
//		{
//			int32_t z = layerNr * thickness + initial;
//			if (z < minZ)
//				continue;
//			if (layerNr < 0)
//				continue;
//
//			SlicerSegment s0;
//			Point3 a0 = p0;
//			Point3 a1 = p1;
//			Point3 a2 = p3;
//
//			bool haveS0 = true;
//			if (a0.z < z && a1.z >= z && a2.z >= z)
//				s0 = project2D(a0, a2, a1, z);
//			else if (a0.z > z && a1.z < z && a2.z < z)
//				s0 = project2D(a0, a1, a2, z);
//
//			else if (a1.z < z && a0.z >= z && a2.z >= z)
//				s0 = project2D(a1, a0, a2, z);
//			else if (a1.z > z && a0.z < z && a2.z < z)
//				s0 = project2D(a1, a2, a0, z);
//
//			else if (a2.z < z && a1.z >= z && a0.z >= z)
//				s0 = project2D(a2, a1, a0, z);
//			else if (a2.z > z && a1.z < z && a0.z < z)
//				s0 = project2D(a2, a0, a1, z);
//			else
//			{
//				//	//Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
//				//	//  on the slice would create two segments			
//				haveS0 = false;
//			}
//
//			SlicerSegment s1;
//			bool haveS1 = true;
//			a0 = p1;
//			a1 = p2;
//			a2 = p3;
//
//			if (a0.z < z && a1.z >= z && a2.z >= z)
//				s1 = project2D(a0, a2, a1, z);
//			else if (a0.z > z && a1.z < z && a2.z < z)
//				s1 = project2D(a0, a1, a2, z);
//
//			else if (a1.z < z && a0.z >= z && a2.z >= z)
//				s1 = project2D(a1, a0, a2, z);
//			else if (a1.z > z && a0.z < z && a2.z < z)
//				s1 = project2D(a1, a2, a0, z);
//
//			else if (a2.z < z && a1.z >= z && a0.z >= z)
//				s1 = project2D(a2, a1, a0, z);
//			else if (a2.z > z && a1.z < z && a0.z < z)
//				s1 = project2D(a2, a0, a1, z);
//			else
//			{
//				//Not all cases create a segment, because a point of a face could create just a dot, and two touching faces
//				//  on the slice would create two segments				
//				haveS1 = false;
//			}
//
//			SlicerSegment s;
//			if (haveS0 && haveS1)
//			{
//				Point diff = s0.end - s1.start;
//				if (shorterThen(diff, 10))
//				{
//					s.start = s0.start;
//					s.end = s1.end;
//					s.faceIndex = faceB;
//				}
//				diff = s1.end - s0.start;
//				if (shorterThen(diff, 10))
//				{
//					s.start = s1.start;
//					s.end = s0.end;
//					s.faceIndex = faceA;
//				}
//			}
//			else if (haveS0)
//			{
//				s = s0;
//				s.faceIndex = faceA;
//			}
//			else if (haveS1)
//			{
//				s = s1;
//				s.faceIndex = faceB;
//			}
//			else
//			{
//				continue;
//			}
//
//			layers[layerNr].faceToSegmentIndex[s.faceIndex] = layers[layerNr].segmentList.size();			
//			s.faceIndex1 = (s.faceIndex == faceA) ? faceB : faceA;
//			s.addedToPolygon = false;
//			layers[layerNr].segmentList.push_back(s);
//		}
//	}
//
//	for (unsigned int layerNr = 0; layerNr<layers.size(); layerNr++)
//	{
//		layers[layerNr].makePolygons(ov, keepNoneClosed, extensiveStitching);
//		//if (layerNr > 0)
//		//	layers[layerNr].MatchingStartPoint(layers[layerNr - 1].GetStartPoint());		
//		//layers[layerNr].Optimize();
//	}
//}

void Slicer::dumpSegmentsToHTML(const char* filename)
{
    float scale = std::max(modelSize.x, modelSize.y) / 1500;
    FILE* f = fopen(filename, "w");
    fprintf(f, "<!DOCTYPE html><html><body>\n");
    for(unsigned int i=0; i<layers.size(); i++)
    {
        fprintf(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" style='width:%ipx;height:%ipx'>\n", int(modelSize.x / scale), int(modelSize.y / scale));
        fprintf(f, "<marker id='MidMarker' viewBox='0 0 10 10' refX='5' refY='5' markerUnits='strokeWidth' markerWidth='10' markerHeight='10' stroke='lightblue' stroke-width='2' fill='none' orient='auto'>");
        fprintf(f, "<path d='M 0 0 L 10 5 M 0 10 L 10 5'/>");
        fprintf(f, "</marker>");
        fprintf(f, "<g fill-rule='evenodd' style=\"fill: gray; stroke:black;stroke-width:1\">\n");
        fprintf(f, "<path marker-mid='url(#MidMarker)' d=\"");
        for(unsigned int j=0; j<layers[i].polygonList.size(); j++)
        {
            PolygonRef p = layers[i].polygonList[j];
            for(unsigned int n=0; n<p.size(); n++)
            {
                if (n == 0)
                    fprintf(f, "M");
                else
                    fprintf(f, "L");
                fprintf(f, "%f,%f ", float(p[n].X - modelMin.x)/scale, float(p[n].Y - modelMin.y)/scale);
            }
            fprintf(f, "Z\n");
        }
        fprintf(f, "\"/>");
        fprintf(f, "</g>\n");
        for(unsigned int j=0; j<layers[i].openPolygonList.size(); j++)
        {
            PolygonRef p = layers[i].openPolygonList[j];
            if (p.size() < 1) continue;
            fprintf(f, "<polyline marker-mid='url(#MidMarker)' points=\"");
            for(unsigned int n=0; n<p.size(); n++)
            {
                fprintf(f, "%f,%f ", float(p[n].X - modelMin.x)/scale, float(p[n].Y - modelMin.y)/scale);
            }
            fprintf(f, "\" style=\"fill: none; stroke:red;stroke-width:1\" />\n");
        }
        fprintf(f, "</svg>\n");
    }
    fprintf(f, "</body></html>");
    fclose(f);
}

void Slicer::Slice(SliceDataStorage *pStorage, HyOutlines *pOutLines)
{
	dumpLayerparts(*pStorage, pOutLines);
	//dumpLayerparts(*pStorage, "c:/slice_test/output_layer.html");
}

SliceDataStorage *Slicer::Slice(SliceDataStorage *pStorage, SModel *pModel, float sliceZ, bool keepNoneClosed, bool extensiveStitching)
{
    SimpleModel* model = NULL;
	if (NULL == pStorage)
	{
		pStorage = new SliceDataStorage();

		model = loadModelFromBuffer(pModel->pBuff, pModel->size);
		if (!model)
		{
			logError("Failed to load from memory.");
			return NULL;
		}

		OptimizedModel* optimizedModel = new OptimizedModel(model, Point3(pModel->posX * 1000, pModel->posY * 1000, -pModel->posZ * 1000));
		//for (unsigned int v = 0; v < model->volumes.size(); v++)
		//{
		//	log("  Face counts: %i -> %i %0.1f%%\n", (int)model->volumes[v].faces.size(), (int)optimizedModel->volumes[v].faces.size(), float(optimizedModel->volumes[v].faces.size()) / float(model->volumes[v].faces.size()) * 100);
		//	log("  Vertex counts: %i -> %i %0.1f%%\n", (int)model->volumes[v].faces.size() * 3, (int)optimizedModel->volumes[v].points.size(), float(optimizedModel->volumes[v].points.size()) / float(model->volumes[v].faces.size() * 3) * 100);
		//}
		delete model;

		//log("Slicing model...\n");
		vector<Slicer*> slicerList;
		for (unsigned int volumeIdx = 0; volumeIdx < optimizedModel->volumes.size(); volumeIdx++)
		{
			//Slicer* slicer = new Slicer(&optimizedModel->volumes[volumeIdx], config.initialLayerThickness - config.layerThickness / 2, config.layerThickness, config.fixHorrible & FIX_HORRIBLE_KEEP_NONE_CLOSED, config.fixHorrible & FIX_HORRIBLE_EXTENSIVE_STITCHING);//org it had a bug related with initail layer thickness
			Slicer* slicer = new Slicer(&optimizedModel->volumes[volumeIdx], pModel->initialLayerThickness * 1000, pModel->layerThickness * 1000, pModel->fixHorrible & FIX_HORRIBLE_KEEP_NONE_CLOSED, pModel->fixHorrible & FIX_HORRIBLE_EXTENSIVE_STITCHING);
			slicerList.push_back(slicer);
			//BHC_DEBUG
			//slicer->dumpSegmentsToHTML("c:/slice_test/output_segment.html");
			for (unsigned int layerNr = 0; layerNr<slicer->layers.size(); layerNr++)
			{
				//Reporting the outline here slows down the engine quite a bit, so only do so when debugging.
				//sendPolygonsToGui("outline", layerNr, slicer->layers[layerNr].z, slicer->layers[layerNr].polygonList);
				//sendPolygonsToGui("openoutline", layerNr, slicer->layers[layerNr].z, slicer->layers[layerNr].openPolygonList);
			}
		}

		pStorage->modelSize = optimizedModel->modelSize;
		pStorage->modelMin = optimizedModel->vMin;
		pStorage->modelMax = optimizedModel->vMax;
		Polygon box;
		box.add(Point(pStorage->modelMin.x, pStorage->modelMax.y));
		box.add(Point(pStorage->modelMin.x, pStorage->modelMin.y));
		box.add(Point(pStorage->modelMax.x, pStorage->modelMin.y));
		box.add(Point(pStorage->modelMax.x, pStorage->modelMax.y));
		pStorage->boundaryBox.add(box);
		//generateSupportGrid(storage.support, optimizedModel, config.supportAngle, config.supportEverywhere > 0, config.supportXYDistance, config.supportZDistance, getRaftThickness());
		//generateSupportGrid(storage.support, optimizedModel, 60, true, 700, 150, 0);
		delete optimizedModel;

		for (unsigned int volumeIdx = 0; volumeIdx < slicerList.size(); volumeIdx++)
		{
			pStorage->volumes.push_back(SliceVolumeStorage());
			createLayerParts(pStorage->volumes[volumeIdx], slicerList[volumeIdx], pModel->fixHorrible & (FIX_HORRIBLE_UNION_ALL_TYPE_A | FIX_HORRIBLE_UNION_ALL_TYPE_B | FIX_HORRIBLE_UNION_ALL_TYPE_C));
			delete slicerList[volumeIdx];
		}
	}
	pModel->pOutLines = new HyOutlines();
	dumpLayerparts(*pStorage, pModel->pOutLines);

	//generateMultipleVolumesOverlap(storage.volumes, 150);
	
	// BHC debug for layer image.
	//dumpLayerparts(*pStorage, "c:/slice_test/output_layer.html");

	// TODO: implement slice by sliceZ value. 
	//for (unsigned int volumeIdx = 0; volumeIdx < storage.volumes.size(); volumeIdx++)
	//{
	//	for (unsigned int layerNr = 0; layerNr<storage.volumes[volumeIdx].layers.size(); layerNr++)
	//	{
	//		int32_t z = 200 + layerNr * 200;//config.initialLayerThickness + layerNr * config.layerThickness;
	//		z += 0;//config.raftBaseThickness + config.raftInterfaceThickness;
	//		SupportPolyGenerator supportGenerator(storage.support, z);
	//	}
	//}



	return pStorage;
}