/** Copyright (C) 2013 Hyvision - Released under terms of the AGPLv3 License */
#include "util.h"
#include <stdio.h>
#include "sliceDataStorage.h"
Vector3 calcNormal(Vector3 p0, Vector3 p1, Vector3 p2)
{
	Vector3 d1 = p0 - p1;
	Vector3 d2 = p1 - p2;
	Vector3 normal = d1.CrossProduct(d2);
	normal.NormalizeSafe();
	return normal;
}

int sign(float value)
{
	if (value < 0)
		return -1;
	else if (value>0)
		return 1;
	return 0;
}

double calcAngle(Vector3 v0, Vector3 v1)
{
	double dotP = (double)v0.DotProduct(v1);
	double sv0 = v0.Length();
	double sv1 = v1.Length();
	double arg = dotP / (sv0 * sv1);
	if (arg > 1.0f)
		arg = 1.0f;
	if (arg < -1.0f)
		arg = -1.0f;

	double angle = acos(arg);
	return RADTODEG(angle);
}

double calcSignedAngle(Vector3 v0, Vector3 v1, Vector3 n)
{	
	double angle = calcAngle(v0, v1);
	double s = sign(n.DotProduct(v0.CrossProduct(v1)));
	// angle in [-179,180]
	double signed_angle = angle * s;

	// angle in [0,360]
	double angle360 = signed_angle;

	return angle360;
}

int compare_PointY(const void* a, const void* b)
{
	int64_t n = (*(Point*)a).Y - (*(Point*)b).Y;
	if (n < 0) return -1;
	if (n > 0) return 1;
	return 0;
}

int compare_CutListData(const void* a, const void* b)
{
	CutListData *pA = (CutListData*)a;
	CutListData *pB = (CutListData*)b;
	int n = pA->pCutList->size() - pB->pCutList->size();
	if (n < 0) return -1;
	if (n > 0) return 1;
	return 0;
}

void appendVector(vector<Point> &dest, vector<Point> src)
{
	for (int i = 0; i < src.size(); i++)
		dest.push_back(src[i]);
}

void findMinMaxY(vector<Point> &points, int &min, int &max)
{
	qsort(points.data(), points.size(), sizeof(Point), compare_PointY);
	if (points.size()>1)
	{
		min = points[0].Y;
		max = points[points.size() - 1].Y;
	} else
	{
		min = -1;
		max = -1;
	}
}

FPoint2 calcDirVector(Point a, Point b)
{
	Point d = b - a;
	FPoint2 dir((double)d.X, (double)d.Y);
	dir = dir.normalize();
	return dir;
}

double calcVectorAngle(Point v0, Point v1)
{
	double dotP = (double)dot(v0, v1);
	double sv0 = (double) vSize(v0);
	double sv1 = (double) vSize(v1);
	double value = (dotP==0|| sv0==0 || sv1==0)?0:dotP / (sv0*sv1);
	if (value < -1.0f)
		value = -1.0f;
	else if (value>1.0f)
		value = 1.0f;
	

	double angle = acos(value);// dotP / (sv0*sv1));
	return RADTODEG(((double)angle));
}

Polygon makeLine(Point p, FPoint2 dir, int dashSize)
{
	Point end(p.X + dir.x*dashSize, p.Y + dir.y*dashSize);
	Polygon polygon;
	polygon.add(p);
	polygon.add(end);
	return polygon;
}

Polygons makeDashedPattern(Point p0, Point p1, int stride, int dashSize, int &lastLen)
{
	if (stride==0)
		stride = dashSize*2;
	Polygons out;			
	Point d = p1 - p0;	
	FPoint2 dp(double(d.X), double(d.Y));
	double size = dp.vSize();
	dp = dp.normalize();

	if (lastLen+size>stride)
	{
		int startLen = stride - lastLen;
		int nline = (int)(lastLen + size)/stride;
		int rest = (lastLen + size) - stride*nline;
		lastLen = rest;		
		Point a(p0.X+dp.x* startLen, p0.Y+dp.y*startLen);
		for (int i=0;i<nline;i++)
		{
			Polygon line = makeLine(a, dp, dashSize);
			out.add(line);
			a.X = a.X + dp.x*stride;
			a.Y = a.Y + dp.y*stride;
		}		
	} else
	{
		lastLen += size;
	}

	return out;
}

Polygons makePolygonDashedPattern(Polygons polygons, int dashStride, int dashSize)
{
	Polygons pieceOfLines;	
	for (unsigned int i=0;i<polygons.size();i++)
	{		
		PolygonRef points = polygons[i];
		unsigned int size = points.size();
		//put a line at the start point
		if (size>1)
		{			
			FPoint2 dir = calcDirVector(points[0], points[1]);
			Polygon line = makeLine(points[0], dir, dashSize);
			pieceOfLines.add(line);
		}

		//fills dash lines in middle of the line
		int rest=0;		
		for (unsigned int j=0;j<size-1;j++)
		{
			Point p0 = points[j];			
			if (j+1<size)
			{
				Point p1 = points[j+1];
				Polygons dashedLines = makeDashedPattern(p0, p1, dashStride, dashSize, rest);
				if (dashedLines.size())
					pieceOfLines.add(dashedLines);
			}
		}

		//automatic closing path
		if (size>2)
		{
			int32_t dist = vSize(points[size-1] - points[0]);
			Point endP = points[size-1];
			//float angle = calcVectorAngle(points[1]-points[0], points[0]-endP);
			if (dist+rest>dashStride)
			{
				Polygons dashedLines = makeDashedPattern(endP, points[0], dashStride, dashSize, rest);
				if (dashedLines.size())
					pieceOfLines.add(dashedLines);
			}
		}
	}

	return pieceOfLines;
}

//line equation
//y = ax + b
//out_a = a (slope)
//out_b = b (constant)
//
inline LineEqResult CalcLineEquation(Point p1, Point p2, double &out_a, double &out_b)
{
	LineEqResult result(LINEQ_NOEXCEPTION, Point(0, 0));

	double dy = (p2.Y - p1.Y);//if dy==0, then x infinite and y == p1.Y == p2.Y
	double dx = (p2.X - p1.X);//if dx==0, then y infinite and x == p1.X == p2.X
	if (dx == 0 && dy == 0)
	{
		result.status = LINEQ_DOT;
	} else
	if (dx==0)
	{
		result.status = LINEQ_INFINITE_Y;
		result.point.X = p1.X;
	} else
	if (dy==0)
	{
		result.status = LINEQ_INFINITE_X;
		result.point.Y = p1.Y;
	} else
	{
		out_a = (dy / dx);
		out_b = -out_a*p1.X + p1.Y;
	}
	return result;
}

bool finiteLineCollusion(Point p1, Point p2, Point q1, Point q2, Point *pOutCrossPoint)
{
	double m=0, m1=0;
	double b=0, b1=0;
	
	LineEqResult rst1 = CalcLineEquation(p1, p2, m, b);
	LineEqResult rst2 = CalcLineEquation(q1, q2, m1, b1);

	Point crossP;
	//exceptions
	if (rst1.status == LINEQ_DOT || rst2.status == LINEQ_DOT)
	{
		return false;
	} else if (rst1.status == LINEQ_INFINITE_X && rst2.status == LINEQ_INFINITE_X)
	{		
		crossP.X = p1.X;//FIXME: just try it
		crossP.Y = p1.Y;		
	} else if (rst1.status == LINEQ_INFINITE_Y && rst2.status == LINEQ_INFINITE_Y)
	{		
		crossP.X = p1.X;
		crossP.Y = p1.Y;//FIXME: just try it		
	} else if (rst1.status== LINEQ_INFINITE_X && rst2.status == LINEQ_NOEXCEPTION)
	{		
		crossP.X = (p1.Y-b1)/m1;
		crossP.Y = p1.Y;		
	} else if (rst1.status == LINEQ_INFINITE_Y && rst2.status == LINEQ_NOEXCEPTION)
	{		
		crossP.X = p1.X;
		crossP.Y = m1*p1.X + b1;		
	} else if (rst2.status == LINEQ_INFINITE_X && rst1.status == LINEQ_NOEXCEPTION)
	{		
		crossP.X = (q1.Y - b) / m;
		crossP.Y = q1.Y;		
	} else if (rst2.status == LINEQ_INFINITE_Y && rst1.status == LINEQ_NOEXCEPTION)
	{		
		crossP.X = q1.X;
		crossP.Y = m*q1.X + b;		
	} else if (rst1.status == LINEQ_INFINITE_X && rst2.status == LINEQ_INFINITE_Y)
	{		
		crossP.X = q1.X;
		crossP.Y = p1.Y;		
	} else if (rst1.status == LINEQ_INFINITE_Y && rst2.status == LINEQ_INFINITE_X)
	{		
		crossP.X = p1.X;
		crossP.Y = q1.Y;		
	} 
	else
	{
		double db = (b1 - b);
		double dm = (m - m1);
		double x = db / dm;
		double y = m*x + b;
		crossP.X = (int)x;
		crossP.Y = (int)y;		
	}

	//text if it is in the range.
	Point max, min;
	max.X = (p1.X > p2.X) ? p1.X : p2.X;
	max.Y = (p1.Y > p2.Y) ? p1.Y : p2.Y;
	min.X = (p1.X < p2.X) ? p1.X : p2.X;
	min.Y = (p1.Y < p2.Y) ? p1.Y : p2.Y;
	if (crossP.X > max.X || crossP.X < min.X || crossP.Y>max.Y || crossP.Y<min.Y)
		return false;

	max.X = (q1.X > q2.X) ? q1.X : q2.X;
	max.Y = (q1.Y > q2.Y) ? q1.Y : q2.Y;
	min.X = (q1.X < q2.X) ? q1.X : q2.X;
	min.Y = (q1.Y < q2.Y) ? q1.Y : q2.Y;
	if (crossP.X > max.X || crossP.X < min.X || crossP.Y>max.Y || crossP.Y < min.Y)
		return false;

	if (pOutCrossPoint)
	{
		//pOutCrossPoint->X = int(x);
		//pOutCrossPoint->Y = int(y);
		*pOutCrossPoint = crossP;
	}
	return true;
}

bool isContainSamePoint(Polygon &polygon, Point p)
{
	for (int i=0;i<int(polygon.size());i++)
	{
		Point point = polygon[i];
		if (point.X == p.X && point.Y==p.Y)
			return true;
	}
	return false;
}

std::vector<CrossResult> findCorssPointsOnPolygon(PolygonRef rPolygon, Point a, Point b, bool findOne)
{
	std::vector<CrossResult> rst;		
	int nPoints = rPolygon.size();
	if (nPoints<2)
		return rst;

	Point prevP = rPolygon[rPolygon.size() - 1];
	for (int j = 0; j<nPoints; j++)
	{
		//it's opened line. it contains 2 points only because of a line. 			
		Point p = rPolygon[j];
		Point crossPoint;
		if (finiteLineCollusion(a, b, prevP, p, &crossPoint))
		{
			rst.push_back(CrossResult(crossPoint, j));
			if (findOne)
				return rst;
		}
		prevP = p;
	}	
	return rst;
}

int findCrossPoints(Polygons &polygons, Point a, Point b, Polygon *pOut, bool instantReturnIfFind)
{	
	std::vector<CrossResult> rst;
	for (int i=0;i<int(polygons.size());i++)
	{
		std::vector<CrossResult> out = findCorssPointsOnPolygon(polygons[i], a, b, instantReturnIfFind);
		rst.insert(rst.end(), out.begin(), out.end());
		if (instantReturnIfFind)
		{
			if (pOut)
			{
				for (int k = 0; k < rst.size(); k++)
					pOut->add(rst[k].p);
			}
			return rst.size();
		}
	}

	if (rst.size() > 0 && pOut)
	{
		for (int i = 0; i < rst.size(); i++)
			pOut->add(rst[i].p);
	}
	return rst.size();
}

std::vector<CrossResult> findCrossPointsHorizontalOnPolygon(PolygonRef points, int y)
{
	std::vector<CrossResult> rst;
	int minY = 0;
	int maxY = 0;
		
	int nPoints = points.size();
	//it's opened line. it contains 2 points only because of a line. 
	if (nPoints <= 2)
		return rst;
	Point p0 = points[nPoints - 1];
	for (int j = 0; j<nPoints; j++)
	{
		Point p1 = points[j];
		minY = p0.Y;
		maxY = p1.Y;
		if (p0.Y>p1.Y)
		{
			minY = p1.Y;
			maxY = p0.Y;
		}
		//equation of strait line given p0, p1
		//y-y1 = a(x - x1)
		//x = (y-y1)/a+x1
		if (y <= maxY && y >= minY)
		{
			int divisor = (p1.Y - p0.Y);
			int x = (divisor == 0) ? 0 : ((y - p0.Y)*(p1.X - p0.X) / divisor + p0.X);
			rst.push_back(CrossResult(Point(x, y), j));
		}
		p0 = p1;
	}	
	return rst;
}

void findCrossPointsHorizontal(Polygon &out, Polygons &polygons, int y)
{
	for (int i = 0; i<int(polygons.size()); i++)
	{
		std::vector<CrossResult> rst = findCrossPointsHorizontalOnPolygon(polygons[i], y);
		for (int k = 0; k < rst.size(); k++)
			out.add(rst[k].p);
	}
}

std::vector<CrossResult> findCrossPointsVerticalOnPolygon(PolygonRef points, int x)
{
	std::vector<CrossResult> rst;
	int minX = 0;
	int maxX = 0;	
	
	int nPoints = points.size();
	//it's opened line. it contains 2 points only because it is a line.
	if (nPoints <= 2)
		return rst;
	Point p0 = points[nPoints - 1];
	for (int j = 0; j<nPoints; j++)
	{
		Point p1 = points[j];
		minX = p0.X;
		maxX = p1.X;
		if (p0.X>p1.X)
		{
			minX = p1.X;
			maxX = p0.X;
		}
		if (x >= minX && x <= maxX)
		{
			int divisor = (p1.X - p0.X);
			int y = (divisor == 0) ? 0 : p0.Y + (p1.Y - p0.Y) * (x - p0.X) / divisor;
			rst.push_back(CrossResult(Point(x, y), j));
		}
		p0 = p1;
	}	
	return rst;
}

void findCrossPointsVertical(Polygon &out, Polygons &polygons, int x)
{		
	for (int i = 0; i<int(polygons.size()); i++)
	{
		std::vector<CrossResult> rst = findCrossPointsVerticalOnPolygon(polygons[i], x);
		for (int k = 0; k < rst.size(); k++)
			out.add(rst[k].p);
	}	
}

int comparePointX(const void* a, const void* b)
{
    int n = ((ClipperLib::IntPoint*)a)->X - ((ClipperLib::IntPoint*)b)->X;
    if (n < 0) return -1;
    if (n > 0) return 1;
    return 0;
}

bool isInsideHorizontal(PolygonRef polygon, Point p)
{
	if (polygon.size()<2)
		return false;
	qsort(polygon.data(), polygon.size(), sizeof(ClipperLib::IntPoint), comparePointX);
	for(int i=0;i<int(polygon.size())/2;i++)
	{
		Point p0 = polygon[i*2];
		Point p1 = polygon[i*2+1];
		int minX=p0.X;
		int maxX=p1.X;
		if (minX>maxX)
		{
			maxX = p0.X;
			minX = p1.X;
		}		
		if (minX<=p.X && p.X<=maxX)
			return true;
	}
	return false;
}

bool isPointInsideInPolygon(Polygons *pOutline, Point p)
{
	Polygon crossPoints;
	findCrossPointsHorizontal(crossPoints, *pOutline, p.Y);
	return isInsideHorizontal(crossPoints, p);	
}

Point findClosestPoint(PolygonRef polygon, Point p)
{
	int minDist = INT_MAX;
	Point minP(-1, -1);
	if (polygon.size()==0)
		return p;
	for(int i=0;i<int(polygon.size());i++)
	{
		Point a = polygon[i];
		Point diff = a - p;
		int dist = (vSize(diff));
		if (minDist>dist)
		{ 
			minDist = dist;
			minP = a;
		}
	}
	return minP;
}

Point findClosestPoint(Polygons &polygons, Point p)
{
	Polygon closestList;
	for (int i = 0; i < polygons.size(); i++)
	{
		if (polygons[i].size()>0)
			closestList.add(findClosestPoint(polygons[i], p));
	}
	return findClosestPoint(closestList, p);
}

Point getInsidePoint(Polygon *pCrossPoints, Point p, Point prevP)
{	
	Point rst = p;	
	if (!pCrossPoints || pCrossPoints->size()<2)
		return p;
	if (!isInsideHorizontal(*pCrossPoints, p))
		rst = findClosestPoint(*pCrossPoints, prevP);
	return rst;
}

Point getAdvancePoint(Polygons outline, Point p, Point prevP)
{
	Point rst;		
	Polygon crossPoints;
	findCrossPoints(outline, prevP, p, &crossPoints, false);
	if (crossPoints.size()==0)
		return p;	
	return findClosestPoint(crossPoints, prevP);
}

void makePolygon(Polygon &out, int64_t x, std::vector<Point> &yarray)
{
	for (int i = 0; i<yarray.size(); i++)
		out.add(Point(x, yarray[i].Y));
}

void optimizePolyline(Polygons &optimizedPolygons, Polygons &polygons)
{	
	PolygonRef optpoly = optimizedPolygons.newPoly();
	for (int i = 0; i < polygons.size(); i++)
	{
		PolygonRef polygon = polygons[i];
		Point prevPoint = polygon[0];
		Point prevDp;
		float prevSlope = -1;

		//optpoly.add(prevPoint);
		for (int j = 1; j < polygon.size()+1; j++)
		{
			int idx = j%polygon.size();
			if (prevPoint.X == polygon[idx].X || prevPoint.Y == polygon[idx].Y)
				continue;

			bool slopeTransition = false;
			Point dp = polygon[idx] - prevPoint;
			float slope = 0;
			if (dp.X == 0 || dp.Y==0)
			{
				slope = 0;
				if (prevSlope == 0 && ((dp.X == 0 && prevDp.Y == 0) || (dp.Y == 0 && prevDp.X == 0)))				
					slopeTransition = true;					
			} else
			{
				slope = (float)dp.Y / (float)dp.X;
				if (slope != prevSlope)
					slopeTransition = true;				
			}		

			if (slopeTransition)
			{				
				optpoly.add(prevPoint);//new point
			} else
			{
				//merge with previouse point
			}			
			//현재 
			//1. 기울기를 구하고 이전값과 같으면이 부분은 확장이다.
			//이전값이 첫번째가 아니면 지운다.
			//2. 이전값과 다르면 이부분은 새로운 시작이다. 새로운 시작 인덱스를 어딘가에 표시 해놓는다.
			//3. 
			prevSlope = slope;
			prevPoint = polygon[idx];			
			prevDp = dp;
		}
	}
}

struct CutData
{
	int index;
	Point crossPoint;
	CutData(int idx, Point p) :index(idx), crossPoint(p.X, p.Y)
	{
	}
};

void extractInsidePolygon(Polygons &out, Polygons &outline, PolygonRef poly)
{
	//1.make cross point list
	vector<CutData> cutList;
	
	for (int i = 0; i < poly.size()-1; i++)
	{
		Polygon crossPoints;
		Point a = poly[i];
		Point b = poly[i + 1];
		findCrossPoints(outline, a, b, &crossPoints, false);
		for (int k = 0; k < crossPoints.size(); k++)		
			cutList.push_back(CutData(i, crossPoints[k]));		
	}

	//2. make piece of polygon that are in the outline.	
	if (cutList.size() == 0)
		return;
	for (int i = 0; i < cutList.size()-1; i++)
	{
		CutData &begin = cutList[i];//begin
		CutData &end = cutList[i + 1];//end		

		if (begin.index == end.index)
		{
			//int32_t advanceSize = 50;
			//Point dir = begin.crossPoint - poly[begin.index];
			//FPoint2 normal((double)dir.X, (double)dir.Y);
			//normal = normal.normalize();
			//Point advancedP(begin.crossPoint.X + normal.x * advanceSize, begin.crossPoint.Y + normal.y * advanceSize);			
			Point center = begin.crossPoint + end.crossPoint;
			center.X /= 2.0f;
			center.Y /= 2.0f;
			Point advancedP = center;
			if (isPointInsideInPolygon(&outline, advancedP))
			{
				PolygonRef piece = out.newPoly();
				piece.add(begin.crossPoint);
				piece.add(end.crossPoint);
			}
		} else
		if (end.index>begin.index)
		{
			Point advancedP = poly[begin.index + 1];			
			if (begin.index + 1 < poly.size() && isPointInsideInPolygon(&outline, advancedP))
			{				
				PolygonRef piece = out.newPoly();
				piece.add(begin.crossPoint);//add the begun point
				for (int k = begin.index + 1; k <= end.index; k++)
					piece.add(poly[k]);
				piece.add(end.crossPoint);//add the end point
			}
		}		
	}	
}

void generateZigzagPattern1(Polygons &patterns, Polygons &outline, int64_t x, int64_t y0, int64_t y1, 
	int lineSpacing, int extrusionWidth, int verticalStride, PointMatrix &matrix, ZIGZAG_ALIGN align)
{
	int64_t diff = y1 - y0;
	int64_t nLine = diff / verticalStride;
	int64_t yspace = diff / nLine;

	//default aligned center	
	int64_t xLeft = x - lineSpacing*0.5f;
	int64_t xRight = x + lineSpacing*0.5f;
	if (align == ZIGZAGALIGN_LEFT)
	{
		xLeft = x;
		xRight = x + (lineSpacing);
	}
	else
	if (align == ZIGZAGALIGN_RIGHT)
	{
		xLeft = x - (lineSpacing);
		xRight = x;
	}

	Polygon poly;
	int64_t y = y0;
	//start point
	Point begin(x, y0);
	poly.add(begin);
	Point prevP = begin;

	//cstyle figure out y0 cross point that pattern goes out of the outline.
	//it can be controlled easily. but in some cases 
	nLine--;
	y += yspace;
	for (int i = 0; i<nLine; i++)
	{
		int64_t xL = (i % 2 == 0) ? xLeft : xRight;
		int64_t xR = (i % 2 == 0) ? xRight : xLeft;
		Point pt(xL, y);		
		poly.add(pt);
		prevP = pt;
		pt = Point(xR, y);		
		prevP = pt;
		poly.add(pt);
		y += yspace;
	}
	//end point
	poly.add(Point(x, y));

	extractInsidePolygon(patterns, outline, poly);	
	//patterns.add(poly);//debug mode
	patterns.unapplyMatrix(matrix);	
}

void generateZigzagPattern(Polygon &poly, Polygons &outline, int64_t x, int64_t y0, int64_t y1, int lineSpacing, int extrusionWidth, int verticalStride, PointMatrix &matrix, ZIGZAG_ALIGN align)
{
	int64_t diff = y1 - y0;
	int64_t nLine = diff / verticalStride;
	int64_t yspace = diff / nLine;

	//default aligned center	
	int64_t xLeft = x - lineSpacing*0.5f;
	int64_t xRight = x + lineSpacing*0.5f;
	if (align == ZIGZAGALIGN_LEFT)
	{
		xLeft = x;
		xRight = x + (lineSpacing);
	} else
	if (align == ZIGZAGALIGN_RIGHT)
	{
		xLeft = x - (lineSpacing);
		xRight = x;
	}
	int64_t y = y0;
	//start point
	Point begin(x, y0);
	poly.add(begin);
	Point prevP = begin;

	//cstyle figure out y0 cross point that pattern goes out of the outline.
	//it can be controlled easily. but in some cases 
	nLine--;
	y += yspace;
	for (int i = 0; i<nLine; i++)
	{
		int64_t xL = (i % 2 == 0) ? xLeft : xRight;
		int64_t xR = (i % 2 == 0) ? xRight : xLeft;
		Point pt(xL, y);
		Polygon crossPoints;
		findCrossPointsHorizontal(crossPoints, outline, y);
		pt = getInsidePoint(&crossPoints, pt, prevP);
		poly.add(pt);
		prevP = pt;
		pt = Point(xR, y);
		pt = getInsidePoint(&crossPoints, pt, prevP);
		prevP = pt;
		poly.add(pt);
		y += yspace;
	}
	//end point
	poly.add(Point(x, y));
	

	//filtering
	//cstyle filtering
	//front
	vector<int> removelist;
	for (int i = 0; i<poly.size(); i++)
	{
		if (i + 1<poly.size() && !isPointInsideInPolygon(&outline, poly[i]))
		{
			Polygon cross;
			if (findCrossPoints(outline, poly[i], poly[i + 1], &cross, true))
			{
				poly[i] = cross[0];
				break;
			} else
				removelist.push_back(i);
		}
	}
	for (int i = 0; i < removelist.size(); i++)
		poly.remove(0);

	//back
	removelist.clear();
	for (int i = poly.size() - 1; i >= 0; i--)
	{
		if (i - 1 >= 0 && !isPointInsideInPolygon(&outline, poly[i]))
		{
			Polygon cross;
			if (findCrossPoints(outline, poly[i], poly[i - 1], &cross, true))
			{
				poly[i] = cross[0];
				break;
			} else
				removelist.push_back(i);
		}
	}
	for (int i = 0; i < removelist.size(); i++)
		poly.remove(poly.size() - 1);
	
	//restore matrix
	for (int i = 0; i < poly.size(); i++)
		poly[i] = matrix.unapply(poly[i]);
}

void extendAABB(AABB &boundary, int offset)
{
	boundary.min.X -= offset;
	boundary.min.Y -= offset;
	boundary.max.X += offset;
	boundary.max.Y += offset;
}

void makeAABBOutline(PolygonRef outPoly, AABB boundary)
{	
	outPoly.add(Point(boundary.min.X, boundary.max.Y));
	outPoly.add(Point(boundary.min.X, boundary.min.Y));
	outPoly.add(Point(boundary.max.X, boundary.min.Y));
	outPoly.add(Point(boundary.max.X, boundary.max.Y));
}

vector< vector<int64_t> > makeCutlist(AABB &boundary, Polygons &outline, int lineSpacing, int offset)
{
	boundary.min.X = ((boundary.min.X / lineSpacing) - 1) * lineSpacing;
	int lineCount = (boundary.max.X - boundary.min.X + (lineSpacing - 1)) / lineSpacing;
	vector<vector<int64_t> > cutList;
	for (int n = 0; n<lineCount; n++)
		cutList.push_back(vector<int64_t>());

	for (unsigned int polyNr = 0; polyNr < outline.size(); polyNr++)
	{
		PolygonRef partOutline = outline[polyNr];
		if (partOutline.size() < 2)
			continue;
		Point p1 = partOutline[partOutline.size() - 1];
		for (unsigned int i = 0; i < outline[polyNr].size(); i++)
		{
			Point p0 = outline[polyNr][i];
			int idx0 = (p0.X - boundary.min.X) / lineSpacing;
			int idx1 = (p1.X - boundary.min.X) / lineSpacing;
			int64_t xMin = p0.X, xMax = p1.X;
			if (p0.X > p1.X) { xMin = p1.X; xMax = p0.X; }
			if (idx0 > idx1) { int tmp = idx0; idx0 = idx1; idx1 = tmp; }
			for (int idx = idx0; idx <= idx1; idx++)
			{
				int x = (idx * lineSpacing) + boundary.min.X + lineSpacing / 2+offset;
				if (x < xMin) continue;
				if (x >= xMax) continue;
				int y = p0.Y + (p1.Y - p0.Y) * (x - p0.X) / (p1.X - p0.X);
				cutList[idx].push_back(y);
			}
			p1 = p0;
		}
	}
	return cutList;
}

vector< vector<Point> > makeCutlistXY(AABB boundary, Polygons &outline, int lineSpacing, int offset)
{
	boundary.min.X = ((boundary.min.X / lineSpacing) - 1) * lineSpacing;
	int lineCount = (boundary.max.X - boundary.min.X + (lineSpacing - 1)) / lineSpacing;

	vector<vector<Point> > cutList;
	for (int n = 0; n<lineCount; n++)
		cutList.push_back(vector<Point>());

	for (unsigned int polyNr = 0; polyNr < outline.size(); polyNr++)
	{
		PolygonRef partOutline = outline[polyNr];
		if (partOutline.size() < 2)
			continue;
		Point p1 = partOutline[partOutline.size() - 1];
		for (unsigned int i = 0; i < outline[polyNr].size(); i++)
		{
			Point p0 = outline[polyNr][i];
			int idx0 = (p0.X - boundary.min.X) / lineSpacing;
			int idx1 = (p1.X - boundary.min.X) / lineSpacing;
			int64_t xMin = p0.X, xMax = p1.X;
			if (p0.X > p1.X) { xMin = p1.X; xMax = p0.X; }
			if (idx0 > idx1) { int tmp = idx0; idx0 = idx1; idx1 = tmp; }
			for (int idx = idx0; idx <= idx1; idx++)
			{
				int x = (idx * lineSpacing) + boundary.min.X + lineSpacing / 2 + offset;
				if (x < xMin) continue;
				if (x >= xMax) continue;
				int y = p0.Y + (p1.Y - p0.Y) * (x - p0.X) / (p1.X - p0.X);
				cutList[idx].push_back(Point(x, y));
			}
			p1 = p0;
		}
	}
	return cutList;
}

void makeTriCutList(TriCutList &outList, AABB boundary, Polygons &outline, int lineSpacing, int offset)
{	
	outList.leftCutList = makeCutlistXY(boundary, outline, lineSpacing, -(lineSpacing - offset) / 2);
	outList.rightCutList = makeCutlistXY(boundary, outline, lineSpacing, (lineSpacing - offset) / 2);
	outList.centerCutList = makeCutlistXY(boundary, outline, lineSpacing);
	outList.CalcBestCase();
}

struct PathCmpData
{
	PolygonRef refPolygon;
	bool	   done;

	PathCmpData(PolygonRef refPoly) :refPolygon(refPoly), done(false)
	{

	}
};

void optimizePolylinePath(Polygons &optPolygons, Polygons &polygons, Point startPoint)
{
	Point startP = startPoint;
	vector<PathCmpData> tmpList;
	for (int i = 0; i < polygons.size(); i++)
	{
		if (polygons[i].size()>1)		
			tmpList.push_back(PathCmpData(polygons[i]));		
	}

	int bestPartIndex = -1;
	int isFront = true;
	int checkCnt = 0;
	vector<int> optimizedIndices;
	while (checkCnt!=tmpList.size())
	{
		int32_t minDist = INT_MAX;
		for (int i = 0; i < tmpList.size(); i++)
		{
			if (!tmpList[i].done)
			{
				PolygonRef part = tmpList[i].refPolygon;
				//compares with the start point
				int32_t dist = vSize(part[0] - startP);
				if (dist < minDist)
				{
					bestPartIndex = i;
					minDist = dist;
					isFront = true;
				}
				//compares with the end point
				dist = vSize(part[part.size() - 1] - startP);
				if (dist < minDist)
				{
					bestPartIndex = i;
					minDist = dist;
					isFront = false;
				}
			}
		}
		
		if (bestPartIndex != -1)
		{
			optimizedIndices.push_back(bestPartIndex);
			tmpList[bestPartIndex].done = true;
			PolygonRef bestPoly = tmpList[bestPartIndex].refPolygon;
			if (!isFront)
				tmpList[bestPartIndex].refPolygon.reverse();
			startP = bestPoly[bestPoly.size() - 1];
			checkCnt++;
		}
	}

	//finaly we get optimized list. So rearranges polygon parts by optimizedIndices
	for (int i = 0; i < optimizedIndices.size(); i++)	
		optPolygons.add(tmpList[optimizedIndices[i]].refPolygon);
}

Polygons makeHeadCleaningPath(int minX, int maxX, int y, int count)
{
	//clean 
	Polygons headCleaning;	
	Polygon cleaningLine;

	int py = y;
	int incY = 1000;
	for (int i=0;i<count;i++)
	{		
		cleaningLine.add(Point(minX, py));
		cleaningLine.add(Point(maxX, py));
		headCleaning.add(cleaningLine);
		py+=incY;
	}
	return headCleaning;
}

int readjustFanspeed(int fanSpeed, int fanThreshold)
{
	//cstyle. there are many kind of fan motor.  generally it is controlled by pwm motor controller. 
	//but some motor would not work at speed 1 or low speed, because of it's own threshold voltage which enables the motor to rotate.
	//To solve this problem, I applied the fanspeedThreshold in the SettingConfig.	
	float maxSpeed = 255.0;
	float threshold = 100 * float(fanThreshold) / maxSpeed;
	float validRange = 100 - threshold;
	float fanSpeedRatio = float(fanSpeed) / 100;
	fanSpeed = threshold + validRange * fanSpeedRatio;
	
	return fanSpeed;
}

void SetPoint(Point &point, int x, int y)
{
	point.X = x;
	point.Y = y;
}

Vector2 RotateVector(float radian, Vector2 p)
{
	Vector2 rst;
	rst.X = cos(radian)*p.X - sin(radian)*p.Y;
	rst.Y = sin(radian)*p.X + cos(radian)*p.Y;
	return rst;
}

Vector2 FindPerpendicularVector(Point a, Point b, float rad, int distanceInMillimeter)
{
	Point d = b - a;
	Vector2 v((double)d.X, (double)d.Y);
	Vector2 c = RotateVector(rad, v);	
	Vector2 vbp = c.NormalizeSafe()*distanceInMillimeter;
	
	return  Vector2((double)a.X, (double)a.Y) + vbp;//a + cni;// vbp.AsPoint();
}

int FindIndexOfPointFromEnd(PolygonRef polygon, int distance)
{
	int lastIndex = polygon.size() - 1;
	int len = 0;
	Point prevP = polygon[lastIndex];
	int rstIndex = lastIndex;

	//find previous index
	int prevIndex = lastIndex;

	for (int i = lastIndex - 1; i>=0; i--)
	{
		//find reversed inside
		Point p = polygon[i];
		Point d = p - prevP;
		int32_t dist = vSize(d);
		len += dist;
		if (len > distance)
		{
			len -= dist;

			Vector2 vd = normalizeSafe(d);
			int restDist = distance - len;
			Point newP(prevP.X + vd.X*(float)restDist, prevP.Y + vd.Y*(float)restDist);
			polygon.insertAt(prevIndex, newP);
			rstIndex = prevIndex;
			break;
		}
		else if (len == distance)
		{
			rstIndex = i;
			break;
		}

		prevP = p;
		prevIndex = i;
	}
	return rstIndex;
}

int FindIndexOfPointFromBegin(PolygonRef polygon, int distance)
{
	//int lastIndex = polygon.size() - 1;
	int len = 0;
	Point prevP = polygon[0];
	int rstIndex = 0;

	//find previous index
	int prevIndex = 0;

	for (int i = 1; i <polygon.size(); i++)
	{
		//find reversed inside
		Point p = polygon[i];
		Point d = p - prevP;
		int32_t dist = vSize(d);
		len += dist;
		if (len > distance)
		{
			len -= dist;

			Point dv = p - prevP;
			Vector2 d = normalizeSafe(dv);
			Point newP(prevP.X + d.X*(distance - len), prevP.Y + d.Y*(distance - len));
			polygon.insertAt(i, newP);
			rstIndex = i;
			break;
		}
		else if (len == distance)
		{
			rstIndex = i;
			break;
		}

		prevP = p;
		prevIndex = i;
	}
	return rstIndex;
}

void TrimPolygonFromBegin(PolygonRef polygon, int distance)
{
	int index = FindIndexOfPointFromBegin(polygon, distance);
	int popCount = index;
	if (popCount > 0 && popCount < polygon.size())
	{
		for (int i = 0; i < popCount; i++)
			polygon.pop_front();
	}
}

void TrimPolygonFromEnd(PolygonRef polygon, int distance)
{
	int index = FindIndexOfPointFromEnd(polygon, distance);
	int popCount = polygon.size() - 1 - index;
	if (popCount > 0 && popCount < polygon.size())
	{
		for (int i = 0; i < popCount; i++)
			polygon.pop_back();
	}
}

void ReOrderPolygon(PolygonRef polygon, int startIdx)
{
	if (startIdx != 0)
	{
		_Polygon newPolygon;
		for (int i = 0; i < polygon.size(); i++)
			newPolygon.add(polygon[(startIdx + i) % polygon.size()]);
		polygon.set(newPolygon);
	}
}

//This function can be used for closed polygon which contains points more than 2.
void MoveStartPosition(PolygonRef polygon, int offset)
{
	if (polygon.size() <= 2)
		return;

	int lastIdx = polygon.size() - 1;
	//connect end point and start point
	if (polygon[lastIdx].X != polygon[0].X || polygon[lastIdx].Y != polygon[0].Y)	
		polygon.add(polygon[0]);
	if (offset < 0)
	{
		//TrimPolygonFromEnd(polygon, -offset);
		int index = FindIndexOfPointFromEnd(polygon, -offset);
		if (index > 0)
		{
			polygon.remove(0);//remove start point which is alread added to the end of array.
			index--;
			ReOrderPolygon(polygon, index);
		}
	}
	else if (offset > 0)
	{
		//TrimPolygonFromBegin(polygon, offset);
		int index = FindIndexOfPointFromBegin(polygon, offset);
		if (index < polygon.size() - 1)
		{
			polygon.remove(0);//remove start point which is alread added to the end of array.
			index--;
			ReOrderPolygon(polygon, index);
		}
	}
	//polygon.remove(polygon.size() - 1);//to disconned the start and end.
	//polygon.addFront(polygon[polygon.size() - 1]);
}

std::vector<CrossResult> FindPointByAngle(PolygonRef polygon, float angle)
{
	Vector2 v = Vector2(cos(angle), sin(angle));
	v *= 1000000.0f;
	
	Point center = polygon.centerOfMass();
	Point b = center + Point(v.X, v.Y);// v.AsPoint();
	return findCorssPointsOnPolygon(polygon, center, b, false);	
}

int FindBigYIn(std::vector<CrossResult> data)
{
	int rst = -1;
	if (data.size()>0)
	{		
		rst = 0;
		int maxY = data[0].p.Y;
		for (int i = 1; i < data.size(); i++)
		{
			if (data[i].p.Y>maxY)
			{
				rst = i;
				maxY = data[i].p.Y;
			}
		}
	}
	return rst;
}

void RearrangeStartPosition(Polygons &polygons, float angleInDeg)
{
	for (int i = 0; i < polygons.size(); i++)
	{		
		PolygonRef poly = polygons[i];
		int startIndex = -1;
		Point startP;
		if (angleInDeg == 90)
		{
			Point center = poly.centerOfMass();
			std::vector<CrossResult> out = findCrossPointsVerticalOnPolygon(poly, center.X);
			int best = FindBigYIn(out);
			if (best >= 0)
			{
				CrossResult cross = out[best];
				startIndex = cross.index;
				startP = cross.p;
			}
		}
		else
		{
			std::vector<CrossResult> out = FindPointByAngle(poly, DEGTORAD(angleInDeg));
			int best = FindBigYIn(out);
			if (best >= 0)
			{
				CrossResult cross = out[best];
				startIndex = cross.index;
				startP = cross.p;
			}
		}
		
		if (startIndex != -1)
		{
			poly.insertAt(startIndex, startP);
			ReOrderPolygon(poly, startIndex);
			//optimize polygon
			if (poly.size() > 3)
			{
				Point v = poly[poly.size() - 1] - poly[0];
				int32_t dist = vSize(v);
				if (dist < 100)
					poly.pop_back();
			}
		}
	}
}

void RearrangeOutlineStartPosition(SliceVolumeStorage& storage, float angleInDeg)
{
	for (unsigned int layerNr = 0; layerNr < storage.layers.size(); layerNr++)
	{				
		for (int i = 0; i < storage.layers[layerNr].parts.size(); i++)		
			RearrangeStartPosition(storage.layers[layerNr].parts[i].outline, angleInDeg);
	}
}

//manual path line
Polygons generateHomeLine(Polygons &polygons, int32_t maxDistance, Point from)
{
	Point endP = findClosestPoint(polygons, from);
	Polygons homeLine;
	PolygonRef poly = homeLine.newPoly();
	Point a = from;
	Point b = endP;
	int32_t len = vSize(b);
	int32_t dist = maxDistance;
	if (len > dist)
	{
		FPoint2 dir = calcDirVector(b, a);
		a.X = b.X + dir.x * dist;
		a.Y = b.Y + dir.y * dist;
	}
	//cstyle force extrusion.
	a.X = 0;
	a.Y = 0;
	poly.add(a);
	poly.add(b);
	return homeLine;
}

Polygons generateHomeLine1(Polygons &polygons, ConfigSettings& config)
{	
	Point home, a, b, c;	
	Polygons homeLine;
	const int marginX = config.homeline_marginX;
	const int marginY = config.homeline_marginY;
	const int startX = config.homeline_startX;
	const int startY = config.homeline_startY;
	const int bedWidth = config.bed_width;
	const int bedHeight = config.bed_height;
	const int centerX = startX + bedWidth/2;
	const int centerY = startY + bedHeight/2;
	home.X = 0;
	home.Y = 0;

	int type = 2;
	if (type == 0 || type == 1)
	{		
		if (type == 0)		
			SetPoint(a, centerX, 0);		
		else
		if (type == 1)		
			SetPoint(a, 0, centerY);		
		Point endP = findClosestPoint(polygons, a);		
		PolygonRef poly = homeLine.newPoly();				
		poly.add(home);
		poly.add(a);
		poly.add(endP);
	} else
	if (type == 2)
	{		
		SetPoint(a, 0, centerY);		
		SetPoint(b, startX + marginX, centerY);		
		SetPoint(c, startX + marginX, startY+marginY);// +(margin < 0) ? 0 : margin);
		Point endP = findClosestPoint(polygons, c);
		PolygonRef poly = homeLine.newPoly();	
		
		poly.add(home);
		poly.add(a);
		poly.add(b);
		poly.add(c);
		poly.add(endP);
	}
	
	return homeLine;
}

Polygons genInitPath(ConfigSettings& config, int32_t distance, Point startP, Polygons *pPolygons)
{
	if (pPolygons)
		return generateHomeLine1(*pPolygons, config);
	return Polygons();
}

void writePolygon(const char *fileName, PolygonRef polygon, const char *comment)
{
	FILE *f = fopen(fileName, "w");
	if (!f)
		return;
	if (comment)
		fprintf(f, "%s \n", comment);
	fprintf(f, "polygon[%d] = {\n", polygon.size());
	for (int i=0;i<polygon.size();i++)
	{
		Point point = polygon[i];
		fprintf(f, "[%d] = %d , %d \n", i, int(point.X), int(point.Y));
	}
	fprintf(f, "}\n");
	fclose(f);
}

void writeOptions(GCodeExport &gcode, ConfigSettings &config, const char *fileName)
{
	gcode.writeLine("");
	gcode.writeComment("**********option description**********");
	if(fileName)
		gcode.writeComment("file: %s", fileName);
	gcode.writeComment("layerThickness=%d", config.layerThickness);
	gcode.writeComment("initialLayerThickness=%d", config.initialLayerThickness);
	gcode.writeComment("filamentDiameter=%d", config.filamentDiameter);
	gcode.writeComment("filamentFlow=%d", config.filamentFlow);
	gcode.writeComment("extrusionWidth=%d", config.extrusionWidth);
	gcode.writeComment("insetCount=%d", config.insetCount);
	gcode.writeComment("downSkinCount=%d", config.downSkinCount);
	gcode.writeComment("upSkinCount=%d", config.upSkinCount);	
	gcode.writeComment("sparseInfillLineDistance=%d", config.sparseInfillLineDistance);
	gcode.writeComment("infillOverlap=%d", config.infillOverlap);
	gcode.writeComment("skirtDistance=%d", config.skirtDistance);
	gcode.writeComment("skirtLineCount=%d", config.skirtLineCount);
	gcode.writeComment("skirtMinLength=%d", config.skirtMinLength);
	
	//Retraction settings
	gcode.writeComment("retractionAmount=%d", config.retractionAmount);
	gcode.writeComment("retractionAmountExtruderSwitch=%d", config.retractionAmountExtruderSwitch);
	gcode.writeComment("retractionSpeed=%d", config.retractionSpeed);
	gcode.writeComment("retractionMinimalDistance=%d", config.retractionMinimalDistance);
	gcode.writeComment("minimalExtrusionBeforeRetraction=%d", config.minimalExtrusionBeforeRetraction);
	gcode.writeComment("retractionZHop=%d", config.retractionZHop);
	gcode.writeComment("retractRestoreCorrection=%d", config.retractRestoreCorrection);

	gcode.writeComment("enableCombing=%d", config.enableCombing);
	gcode.writeComment("enableOozeShield=%d", config.enableOozeShield);
	//gcode.writeComment("wipeTowerSize=%d", config.wipeTowerSize);
	gcode.writeComment("multiVolumeOverlap=%d", config.multiVolumeOverlap);
	
	gcode.writeComment("initialSpeedupLayers=%d", config.initialSpeedupLayers);
	gcode.writeComment("initialLayerSpeed=%d", config.initialLayerSpeed);
	gcode.writeComment("printSpeed=%d", config.printSpeed);
	gcode.writeComment("infillSpeed=%d", config.infillSpeed);
	gcode.writeComment("inset0Speed=%d", config.inset0Speed);
	gcode.writeComment("insetXSpeed=%d", config.insetXSpeed);
	gcode.writeComment("moveSpeed=%d", config.moveSpeed);
	gcode.writeComment("skirtSpeed=%d", config.skirtSpeed);
	gcode.writeComment("startSpeed=%d", config.startSpeed);
	gcode.writeComment("fanFullOnLayerNr=%d", config.fanFullOnLayerNr);
	//Support material
	gcode.writeComment("supportAngle=%d", config.supportAngle);
	gcode.writeComment("supportEverywhere=%d", config.supportEverywhere);
	gcode.writeComment("supportLineDistance=%d", config.supportLineDistance);
	gcode.writeComment("supportXYDistance=%d", config.supportXYDistance);
	gcode.writeComment("supportZDistance=%d", config.supportZDistance);
	gcode.writeComment("supportExtruder=%d", config.supportExtruder);
	
	//Cool settings
	gcode.writeComment("minimalLayerTime=%d", config.minimalLayerTime);
	gcode.writeComment("minimalFeedrate=%d", config.minimalFeedrate);
	gcode.writeComment("coolHeadLift=%d", config.coolHeadLift);
	gcode.writeComment("fanSpeedMin=%d", config.fanSpeedMin);
	gcode.writeComment("fanSpeedMax=%d", config.fanSpeedMax);
	
	//Raft settings
	gcode.writeComment("raftType=%d", config.raftType);
	gcode.writeComment("raftMargin=%d", config.raftMargin);
	gcode.writeComment("raftLineSpacing=%d", config.raftLineSpacing);
	gcode.writeComment("raftBaseThickness=%d", config.raftBaseThickness);
	gcode.writeComment("raftBaseLinewidth=%d", config.raftBaseLinewidth);
	gcode.writeComment("raftInterfaceThickness=%d", config.raftInterfaceThickness);
	gcode.writeComment("raftInterfaceLinewidth=%d", config.raftInterfaceLinewidth);	
	gcode.writeComment("raftPlaneWidth1=%d", config.raftPlaneWidth1);
	gcode.writeComment("raftPlaneWidth2=%d", config.raftPlaneWidth2);
	gcode.writeComment("raftPlaneThickness1=%d", config.raftPlaneThickness1);
	gcode.writeComment("raftPlaneThickness2=%d", config.raftPlaneThickness2);
	gcode.writeComment("raftBaseSpeed=%d", config.raftBaseSpeed);
	gcode.writeComment("raftInterfaceSpeed=%d", config.raftInterfaceSpeed);
	gcode.writeComment("raftAirGap=%d", config.airGap);

	
	gcode.writeComment("matrix=%.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f", config.matrix.m[0][0], config.matrix.m[0][1], config.matrix.m[0][2], 
		config.matrix.m[1][0], config.matrix.m[1][1], config.matrix.m[1][2], 
		config.matrix.m[2][0], config.matrix.m[2][1], config.matrix.m[2][2]);
	gcode.writeComment("objectPosition=%d, %d", config.objectPosition.X, config.objectPosition.Y);
	gcode.writeComment("objectSink=%d", config.objectSink);
	
	gcode.writeComment("fixHorrible=%d", config.fixHorrible);
	gcode.writeComment("spiralizeMode=%d", config.spiralizeMode);
	gcode.writeComment("gcodeFlavor=%d", config.gcodeFlavor);	
	gcode.writeComment("layerFill=%d", config.layerSupport);
	gcode.writeComment("layerSupportTopMarginCount=%d", config.layerSupportTopMarginCount);
	
	//additional option
	//cstyle add dash pattern
	//gcode.writeComment("raftDashStride=%d", config.raftDashStride);
	//gcode.writeComment("raftDashSize=%d", config.raftDashSize);
	//gcode.writeComment("wallInOutOrder=%d", config.wallInOutOrder);	

	//gcode.writeComment("insetEndRetraction=%d", config.insetEndRetraction);
	//wipe head
	//gcode.writeComment("wipetype=%d", config.wipetype);//0, 1
	//gcode.writeComment("wipePosition=%d", config.wipePosition.X);
	//gcode.writeComment("wipePosition1=%d", config.wipePosition1.X);

	//gcode.writeComment("wipeExtrusionAmount=%d", config.wipeExtrusionAmount);

	gcode.writeComment("homeline_marginX=%d", config.homeline_marginX);
	gcode.writeComment("homeline_marginY=%d", config.homeline_marginY);
	gcode.writeComment("homeline_startX=%d", config.homeline_startX);
	gcode.writeComment("homeline_startY=%d", config.homeline_startY);
	gcode.writeComment("bed_width=%d", config.bed_width);
	gcode.writeComment("bed_height=%d", config.bed_height);

	//gcode.writeComment("innerWallEndControlType=%d", config.innerWallEndControlType);
	//gcode.writeComment("innerwallEndReduceRate=%d", config.innerwallEndReduceRate);
	gcode.writeComment("innerWallEndGap=%d", config.innerWallEndGap);
	//gcode.writeComment("outerWallEndGap=%d", config.outerWallEndGap);

	//gcode.writeComment("pullEndOfWallTravelDistance=%d", config.pullEndOfWallTravelDistance);
	//gcode.writeComment("pullEndOfWallAngle=%d", config.pullEndOfWallAngle);
	//gcode.writeComment("travelSpeedMatching=%d", config.travelSpeedMatching);
	//gcode.writeComment("outwallStartOffset=%d", config.outerWallStartOffset);
	gcode.writeComment("fixedWallBeginAngle=%d", config.fixedWallBeginAngle);
	//gcode.writeComment("finalizeExtLift=%d", config.finalizeExtLift);

	gcode.writeComment("finalizeRetraction=%d", config.finalizeRetraction);
	gcode.writeComment("fanspeedThreshold=%d", config.fanspeedThreshold);
	gcode.writeComment("initialLayerFillType=%d", config.initialLayerFillType);
	gcode.writeComment("**********end of the option description**********\n");
}