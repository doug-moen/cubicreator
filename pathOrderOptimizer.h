/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef PATHOPTIMIZER_H
#define PATHOPTIMIZER_H

#include <stdint.h>
#include "utils/polygon.h"

class PathOrderOptimizer
{
public:
    Point startPoint;
    vector<PolygonRef> polygons;
    vector<int> polyStart;
    vector<int> polyOrder;
	//cstyle
	Point closestStartPoint;

    PathOrderOptimizer(Point startPoint)
    {
        this->startPoint = startPoint;
		this->closestStartPoint = startPoint;
    }

    void addPolygon(PolygonRef polygon)
    {
        this->polygons.push_back(polygon);
    }
    
    void addPolygons(Polygons& polygons)
    {
        for(unsigned int i=0;i<polygons.size(); i++)
            this->polygons.push_back(polygons[i]);
    }
    
    void optimize(bool keepPointOrder=false);
	void makePathOrder();
};

#endif//PATHOPTIMIZER_H
