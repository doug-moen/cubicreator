/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef INFILL_H
#define INFILL_H

#include "utils/polygon.h"
#include "settings.h"
enum PATTERN_TYPE
{
	PATTERN_LINE,
	PATTERN_ZIGZAG,
	PATTERN_DASH
};

void generateConcentricInfill(Polygons outline, Polygons& result, int offsets[], int offsetsSize);
void generateLineInfill(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation);
void generateLineInfillEx(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation, 
	PATTERN_TYPE patternType, int arg0 = 0, int arg1 = 0, Polygons *pRefBoundary = NULL, Polygons *pRefInsets = NULL);
#endif//INFILL_H
