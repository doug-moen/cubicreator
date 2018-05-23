/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include "polygonOptimizer.h"
#include "settings.h"
#include "utils/util.h"

void optimizePolygon(PolygonRef poly)
{
    Point p0 = poly[poly.size()-1];
    for(unsigned int i=0;i<poly.size();i++)
    {
        Point p1 = poly[i];
        //if (shorterThen(p0 - p1, 10))
		if (shorterThen(p0 - p1, MINIMAL_POLYGON_LENGTH))
        {
            poly.remove(i);
            i --;
        }
		else
		{
            Point p2;
            if (i < poly.size() - 1)
                p2 = poly[i+1];
            else
                p2 = poly[0];
            
			Point p10 = p1 - p0;
			Point p12 = p1 - p2;

            Point diff0 = normal(p10, 1000000);
            Point diff2 = normal(p12, 1000000);
            int64_t d = dot(diff0, diff2);
            if (d  < -999999000000LL)			
			//double angle = abs(calcVectorAngle(p10, p12));			
			//if (179.99999f<angle && angle<180.00001f)
			//if (179.999999f<angle && angle<180.00001f)
				//double diffAngle = abs(180.0f - angle);
				//if (diffAngle>0.001f)
            {
                poly.remove(i);
                i --;
            }
			else
			{
                p0 = p1;
            }			
        }
    }
}

void optimizePolygons(Polygons& polys)
{
    for(unsigned int n=0;n<polys.size();n++)
    {
        optimizePolygon(polys[n]);
        if (polys[n].size() < 3)
        {
            polys.remove(n);
            n--;
        }
    }
}
