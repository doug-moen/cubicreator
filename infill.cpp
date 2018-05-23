/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include "infill.h"
#include "utils/util.h"

void generateConcentricInfill(Polygons outline, Polygons& result, int offsets[], int offsetsSize)
{
    int step = 0;
    while(1)
    {
        for(unsigned int polygonNr=0; polygonNr<outline.size(); polygonNr++)
            result.add(outline[polygonNr]);
        outline = outline.offset(-offsets[step]);
        if (outline.size() < 1)
            break;
        step = (step + 1) % offsetsSize;
    }
}

int compare_int64_t(const void* a, const void* b)
{
    int64_t n = (*(int64_t*)a) - (*(int64_t*)b);
    if (n < 0) return -1;
    if (n > 0) return 1;
    return 0;
}

//int compare_PointY(const void* a, const void* b)
//{
//	int64_t n = (*(Point*)a).Y - (*(Point*)b).Y;
//	if (n < 0) return -1;
//	if (n > 0) return 1;
//	return 0;
//}

void generateLineInfill(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation)
{
    Polygons outline = in_outline.offset(extrusionWidth * infillOverlap / 100);
    PointMatrix matrix(rotation);
    
    outline.applyMatrix(matrix);
    
    AABB boundary(outline);
    
    boundary.min.X = ((boundary.min.X / lineSpacing) - 1) * lineSpacing;
    int lineCount = (boundary.max.X - boundary.min.X + (lineSpacing - 1)) / lineSpacing;
    vector<vector<int64_t> > cutList;
    for(int n=0; n<lineCount; n++)
        cutList.push_back(vector<int64_t>());

    for(unsigned int polyNr=0; polyNr < outline.size(); polyNr++)
    {
        Point p1 = outline[polyNr][outline[polyNr].size()-1];
        for(unsigned int i=0; i < outline[polyNr].size(); i++)
        {
            Point p0 = outline[polyNr][i];
            int idx0 = (p0.X - boundary.min.X) / lineSpacing;
            int idx1 = (p1.X - boundary.min.X) / lineSpacing;
            int64_t xMin = p0.X, xMax = p1.X;
            if (p0.X > p1.X) { xMin = p1.X; xMax = p0.X; }
            if (idx0 > idx1) { int tmp = idx0; idx0 = idx1; idx1 = tmp; }
            for(int idx = idx0; idx<=idx1; idx++)
            {
                int x = (idx * lineSpacing) + boundary.min.X + lineSpacing / 2;
                if (x < xMin) continue;
                if (x >= xMax) continue;
                int y = p0.Y + (p1.Y - p0.Y) * (x - p0.X) / (p1.X - p0.X);
                cutList[idx].push_back(y);
            }
            p1 = p0;
        }
    }
    
    int idx = 0;
    for(int64_t x = boundary.min.X + lineSpacing / 2; x < boundary.max.X; x += lineSpacing)
    {
        qsort(cutList[idx].data(), cutList[idx].size(), sizeof(int64_t), compare_int64_t);
        for(unsigned int i = 0; i + 1 < cutList[idx].size(); i+=2)
        {
            if (cutList[idx][i+1] - cutList[idx][i] < extrusionWidth / 5)
                continue;
            PolygonRef p = result.newPoly();
            p.add(matrix.unapply(Point(x, cutList[idx][i])));
            p.add(matrix.unapply(Point(x, cutList[idx][i+1])));
        }
        idx += 1;
    }
}

void generateLineInfillEx(const Polygons& in_outline, Polygons& result, int extrusionWidth, int lineSpacing, int infillOverlap, double rotation, 
	PATTERN_TYPE patternType, int arg0, int arg1, Polygons *pRefBoundary, Polygons *pRefInsets)
{
    Polygons outline = in_outline.offset(extrusionWidth * infillOverlap / 100);

    PointMatrix matrix(rotation);
    outline.applyMatrix(matrix);	
	
	AABB boundary;
	if (pRefBoundary)
	{
		Polygons refBoundary = *pRefBoundary;
		refBoundary.applyMatrix(matrix);
		boundary.calculate(refBoundary);
	} else
		boundary.calculate(outline);

	vector<vector<Point> > *pCutList = NULL;

	Polygons testOutline;
	PolygonRef box = testOutline.newPoly();
	extendAABB(boundary, 1200);
	makeAABBOutline(box, boundary);
		
	vector<vector<Point> > cutList;
	if (patternType == PATTERN_ZIGZAG)
	{
		cutList = makeCutlistXY(boundary, testOutline, lineSpacing);
		pCutList = &cutList;		
	} else
	{
		cutList = makeCutlistXY(boundary, outline, lineSpacing);
		pCutList = &cutList;
	}

    int idx = 0;	
	Point prevP(-1, -1);
	int total = pCutList->size();	

	int64_t x;
	bool first = true;
	Polygon p;    
	for (int n = 0; n < total;n++)
    {
		vector<Point> &cutline = (*pCutList)[idx];
		qsort(cutline.data(), cutline.size(), sizeof(Point), compare_PointY);
		for (unsigned int i = 0; i + 1 < cutline.size(); i += 2)
        {
			if (cutline[i + 1].Y - cutline[i].Y < extrusionWidth / 5)
                continue;
			
			x = cutline[i].X;
			int64_t y0 = cutline[i].Y;
			int64_t y1 = cutline[i + 1].Y;

			if (patternType == PATTERN_DASH)
			{
				Polygons lines;
				Polygon line;
				line.add(matrix.unapply(Point(x, y0)));
				line.add(matrix.unapply(Point(x, y1)));
				lines.add(line);
				int dashStride = arg0;
				int dashSize =arg1;
				Polygons dashedPattern = makePolygonDashedPattern(lines, dashStride, dashSize);
				result.add(dashedPattern);
			} else
			if (patternType==PATTERN_ZIGZAG)
			{
				Polygon part;				
				int64_t diff = y1-y0;
				int64_t verticalStride = int64_t(arg0);
				if (diff<verticalStride)
					goto DEFAULT_PATTERN;			
				
				Polygons parts;
				generateZigzagPattern1(parts, outline, x, y0, y1, lineSpacing - extrusionWidth * 2, extrusionWidth, verticalStride, matrix, ZIGZAGALIGN_CENTER);
				result.add(parts);				
			} else
			{
DEFAULT_PATTERN:
				//connect previous line
				if (arg0==1 && prevP.X>-1 && idx && cutList[idx-1].size() && cutList[idx].size() && (idx>0 && cutList[idx].size()<=cutList[idx-1].size()))
				{					
					Point top(x, cutList[idx].front().Y);
					Point bottom(x, cutList[idx].back().Y);
					Point pt = (idx%2==0)?top:bottom;					
					Polygon poly;
					makePolygon(poly, prevP.X, cutList[idx - 1]);
					Point closestP = findClosestPoint(poly, pt);
					//find line collision
					if (pRefInsets)
						pRefInsets->applyMatrix(matrix);
					
					if (pRefInsets && findCrossPoints(*pRefInsets, closestP, pt, NULL, true)==0)
					{
						Polygon bridge;
						bridge.add(matrix.unapply(closestP));
						bridge.add(matrix.unapply(pt));
						result.add(bridge);
					}

					if (arg1 == 1)
					{
						pt = (idx % 2 == 1) ? top : bottom;
						Polygon poly;
						makePolygon(poly, prevP.X, cutList[idx - 1]);
						Point closestP = findClosestPoint(poly, pt);
						//find line collision
						if (pRefInsets)
							pRefInsets->applyMatrix(matrix);

						if (pRefInsets && findCrossPoints(*pRefInsets, closestP, pt, NULL, true) == 0)
						{
							Polygon bridge;
							bridge.add(matrix.unapply(closestP));
							bridge.add(matrix.unapply(pt));
							result.add(bridge);
						}
					}
				}
				Polygon lines;
				lines.add(matrix.unapply(Point(x, y0)));
				lines.add(matrix.unapply(Point(x, y1)));
				result.add(lines);
				prevP.X = x;
				prevP.Y = y1;
			}
        }		
        idx += 1;
    }

	//cstyle test add a additional zigzag pattern line when it has a overed pattern gap at the rest empty space.
	//find last index
	if (patternType == PATTERN_ZIGZAG)
	{
		int endIndex = -1;
		for (int j = cutList.size() - 1; j >= 0; j--)
		{
			if (cutList[j].size())
			{
				endIndex = j;
				break;
			}
		}

		if (endIndex >= 0)
		{
			int subEndIndex = cutList[endIndex].size() - 1;
			if (subEndIndex >= 0)
			{
				Point p = cutList[endIndex][subEndIndex];
				AABB bound(outline);
				int nextX = p.X + lineSpacing / 2 + extrusionWidth;
				if (nextX < bound.max.X)
				{
					Polygon crossPoints;
					findCrossPointsVertical(crossPoints, testOutline, nextX);
					vector<Point> cutPoints;
					for (int i = 0; i < crossPoints.size(); i++)
						cutPoints.push_back(crossPoints[i]);
					qsort(cutPoints.data(), cutPoints.size(), sizeof(Point), compare_PointY);
					
					int64_t verticalStride = int64_t(arg0);
					for (unsigned int i = 0; i + 1 < cutPoints.size(); i += 2)
					{
						if (cutPoints[i + 1].Y - cutPoints[i].Y < extrusionWidth / 5)
							continue;

						x = cutPoints[i].X;
						int64_t y0 = cutPoints[i].Y;
						int64_t y1 = cutPoints[i + 1].Y;
						Polygons parts;
						generateZigzagPattern1(parts, outline, x, y0, y1, lineSpacing - extrusionWidth * 2, extrusionWidth, verticalStride, matrix, ZIGZAGALIGN_LEFT);
						result.add(parts);
					}
				}
			}
		}
	}

}