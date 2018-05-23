#include "GCodeControl.h"
#include "utils\polygon.h"
#include "utils\util.h"

void GCodeControl::ApplyEndGap(PolygonRef polygon)
{
	if (endGap > 0)
	{
		//first make a closed polygon.
		polygon.add(polygon[0]);
		TrimPolygonFromEnd(polygon, endGap);
	}
}

int GCodeControl::ApplyEndReduce(PolygonRef polygon)
{
	int reduceIndex = -1;
	if (endReduce && polygon.size()>1 && polygon.polygonLength()>reduceDist)
	{
		polygon.add(polygon[0]);
		//reduceIndex = FindIndexOfPointFromEnd(polygon, reducePrepareDist);
		reduceIndex = FindIndexOfPointFromEnd(polygon, reduceDist);
		polygon.remove(polygon.size() - 1);
	}
	else
		endReduce = false;
	return reduceIndex;
}