/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include "pathOrderOptimizer.h"

#ifndef DBL_MAX
#define DBL_MAX (0xFFFFFFFFFFFFFFFFLL)
#endif
void PathOrderOptimizer::optimize(bool keepPointOrder)
{
	if (keepPointOrder)
	{
		makePathOrder();
		return;
	}

    std::vector<bool> picked;
	double maxDist = DBL_MAX;	
    for(unsigned int i=0;i<polygons.size(); i++)
    {
        int best = -1;
        float bestDist = maxDist;
        PolygonRef poly = polygons[i];
        for(unsigned int j=0; j<poly.size(); j++)
        {
            float dist = vSize2f(poly[j] - startPoint);
            if (dist < bestDist)
            {
                best = j;
                bestDist = dist;
            }
        }
        polyStart.push_back(best);
        picked.push_back(false);
    }

	//cstyle: 시작점으로 부터 가까운 폴리곤의 연속된 그리기 순서를 찾아낸다.
    Point p0 = startPoint;
    for(unsigned int n=0; n<polygons.size(); n++)
    {
        int best = -1;
        float bestDist = maxDist;
        for(unsigned int i=0;i<polygons.size(); i++)
        {
            if (picked[i] || polygons[i].size() < 1)
                continue;
			size_t polygonSize = polygons[i].size();            
			//below 2 case can be merged.
			if (polygonSize == 2)
            {
                float dist = vSize2f(polygons[i][0] - p0);
                if (dist < bestDist)
                {
                    best = i;
                    bestDist = dist;
                    polyStart[i] = 0;
                }
                dist = vSize2f(polygons[i][1] - p0);
                if (dist < bestDist)
                {
                    best = i;
                    bestDist = dist;
                    polyStart[i] = 1;
                }
            } if (polygonSize>2 && keepPointOrder)//cstyle
			{
				float dist = vSize2f(polygons[i][0] - p0);
                if (dist < bestDist)
                {
                    best = i;
                    bestDist = dist;
                    polyStart[i] = 0;
                }
				unsigned int endIndex = polygonSize-1;
                dist = vSize2f(polygons[i][endIndex] - p0);
                if (dist < bestDist)
                {
                    best = i;
                    bestDist = dist;
                    polyStart[i] = endIndex;
                }
			} else
			{
                float dist = vSize2f(polygons[i][polyStart[i]] - p0);
                if (dist < bestDist)
                {
                    best = i;
                    bestDist = dist;
                }
            }
        }
        if (best > -1)
        {
            if (polygons[best].size() == 2)
            {
                p0 = polygons[best][(polyStart[best] + 1) % 2];//cstyle: it figures out the end point that is the farest from the start point.
            } else//cstyle
			if (polygons[best].size()>2 && keepPointOrder)
			{
				unsigned int endIndex = polygons[best].size()-1;
				if (polyStart[best]==0)
					p0 = polygons[best][endIndex];
				else
					p0 = polygons[best][0];
			} else
			{
                p0 = polygons[best][polyStart[best]];
            }
            picked[best] = true;
            polyOrder.push_back(best);
        }
    }
    
	//cstyle: recalculates polyStart
    p0 = startPoint;
    for (unsigned int n=0; n<polyOrder.size(); n++)
    {
        int nr = polyOrder[n];
        int best = -1;
        float bestDist = maxDist;
		PolygonRef polygon = polygons[nr];
		//cstyle
		if (polygon.size()>2 && keepPointOrder)
		{
			float dist = vSize2f(polygon[0] - p0);
			if (dist < bestDist)
			{
				best = 0;
				bestDist = dist;
			}
			unsigned int endIndex = polygon.size()-1;
			dist = vSize2f(polygon[endIndex] - p0);
			if (dist < bestDist)
			{
				best = endIndex;
				bestDist = dist;
			}
		} else
		{
        for(unsigned int i=0;i<polygons[nr].size(); i++)
			{
				float dist = vSize2f(polygons[nr][i] - p0);
				if (dist < bestDist)
				{
					best = i;
					bestDist = dist;
				}
			}
		}
        polyStart[nr] = best;
        if (polygons[nr].size() <= 2)
        {
            p0 = polygons[nr][(best + 1) % 2];
        } else
		if (polygons[nr].size() > 2 && keepPointOrder)
		{
			unsigned int endIndex = polygons[nr].size()-1;
			if (best==0)
				p0 = polygons[nr][endIndex];
			else
				p0 = polygons[nr][0];
		} else
		{			
            p0 = polygons[nr][best];
        }
    }
	if (polygons.size() && polygons[0].size() && polyOrder.size() && polyStart.size())
		closestStartPoint = polygons[polyOrder[0]][polyStart[polyOrder[0]]];
}

void PathOrderOptimizer::makePathOrder()
{
	for(unsigned int i=0;i<polygons.size(); i++)
    {        
		polyOrder.push_back(i);
		//if (i%2==0)
			polyStart.push_back(0);
		//else
		//	polyStart.push_back(polygons[i]->size()-1);
    }
}
