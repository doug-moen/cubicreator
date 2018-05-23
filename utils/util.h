/** Copyright (C) 2013 Hyvision - Released under terms of the AGPLv3 License */
#ifndef UTIL_H
#define UTIL_H

#include "polygon.h"
#include "floatpoint.h"
#include "gcodeExport.h"
//#include "settings.h"
//#include "AABB.h"
#include <sstream>

//using namespace cura;
#ifdef WIN32
#include <algorithm>
#endif

#define PI			   (3.14159f)
#define DEGTORAD(X)    (((X)*PI)/180)
#define RADTODEG(X)    (((X)*180)/PI)

enum ZIGZAG_ALIGN
{
	ZIGZAGALIGN_CENTER = 0,
	ZIGZAGALIGN_LEFT,
	ZIGZAGALIGN_RIGHT
};

enum LINEQUATION_STATUS
{
	LINEQ_NOEXCEPTION = 0,
	LINEQ_DOT,
	LINEQ_INFINITE_X,
	LINEQ_INFINITE_Y,
};

struct LineEqResult
{
	LINEQUATION_STATUS status;
	Point point;
	LineEqResult(LINEQUATION_STATUS status_, Point p) :status(status_), point(p)
	{
	}
};

struct CrossResult
{
	Point p;
	int index;

	CrossResult(Point _p, int _index)
	{
		p = _p;
		index = _index;
	}
};

class SliceVolumeStorage;

int compare_PointY(const void* a, const void* b);
void appendVector(vector<Point> &dest, vector<Point> src);
void findMinMaxY(vector<Point> &points, int &min, int &max);


struct CutListData
{
	vector< Point>		*pCutList;
	ZIGZAG_ALIGN		myAlign;

	CutListData(vector<Point> *pList, ZIGZAG_ALIGN align) :
	pCutList(pList), myAlign(align)
	{
	}
};

int compare_CutListData(const void* a, const void* b);

struct TriCutList
{	
	vector< vector< Point> >	leftCutList;
	vector< vector< Point> >	rightCutList;
	vector< vector< Point> >	centerCutList;
	vector< vector< Point> >	bestCutList;
	vector< ZIGZAG_ALIGN >		bestAligns;

	TriCutList()
	{
	}

	void CalcBestCase()
	{
		vector<vector<Point> > &cutList = centerCutList;
		int minY = -1, maxY = -1;
		for (int i = 0; i < cutList.size(); i++)
		{			
			bestAligns.push_back(ZIGZAGALIGN_CENTER);
			bestCutList.push_back(vector<Point>());

			vector<Point> &center = centerCutList[i];
			vector<Point> &left = leftCutList[i];
			vector<Point> &right = rightCutList[i];
			if (!(center.size()>1 || left.size() > 1 || right.size() > 1))
				continue;

			//here figure out best case
			vector<vector<Point> *> bests;
			if (center.size()>1 && center.size() == left.size() && left.size() == right.size())
			{				
				bestAligns[i] = ZIGZAGALIGN_CENTER;
				bests.push_back(&center);
				bests.push_back(&left);
				bests.push_back(&right);
			} else			
			{
				vector<CutListData> testList;
				if (center.size()>1)
					testList.push_back(CutListData(&centerCutList[i], ZIGZAGALIGN_CENTER));
				if (left.size()>1)
					testList.push_back(CutListData(&leftCutList[i], ZIGZAGALIGN_LEFT));
				if (right.size()>1)
					testList.push_back(CutListData(&rightCutList[i], ZIGZAGALIGN_RIGHT));
				qsort(testList.data(), testList.size(), sizeof(CutListData), compare_CutListData);				
				bestAligns[i] = testList[0].myAlign;
				for (int k = 0; k < testList.size(); k++)
					bests.push_back(testList[k].pCutList);
			}		

			//베스트를 정렬 해놓는다.
			vector<Point> &bestCut = *bests[0];
			for (int k = 0; k < bests.size();k++)
				qsort(bests[k]->data(), bests[k]->size(), sizeof(Point), compare_PointY);
			
			for (int j = 0; j < bestCut.size() / 2; j++)
			{
				vector<Point> tmpList;
				for (int k = 0; k < bests.size(); k++)
				{
					vector<Point> &cur = (*bests[k]);
					if (bestCut.size() == cur.size())
					{						
						tmpList.push_back(cur[j*2]);
						tmpList.push_back(cur[j*2+1]);
					} else 
					if (j == 0)					
						tmpList.push_back(cur[0]);					
					if (j+1 == bestCut.size() - 1)
						tmpList.push_back(cur[cur.size() - 1]);
					
				}
				int minY, maxY;
				findMinMaxY(tmpList, minY, maxY);			
				bestCut[j * 2].Y = minY;
				bestCut[j*2+1].Y = maxY;				
			}

			vector<Point> myBest = bestCut;
			bestCutList[i] = myBest;
		}
	}
};

Vector3 calcNormal(Vector3 p0, Vector3 p1, Vector3 p2);
int sign(float value);
double calcAngle(Vector3 v0, Vector3 v1);
double calcSignedAngle(Vector3 v0, Vector3 v1, Vector3 n);

FPoint2 calcDirVector(Point a, Point b);
double calcVectorAngle(Point v0, Point v1);
Polygon makeLine(Point p, FPoint2 dir, int dashSize);
Polygons makeDashedPattern(Point p0, Point p1, int stride, int dashSize, int &lastLen);
Polygons makePolygonDashedPattern(Polygons polygons, int dashStride, int dashSize);

bool FiniteLineCollusion(Point p1, Point p2, Point q1, Point q2,Point *pOutCrossPoint);

int findCrossPoints(Polygons &polygons, Point a, Point b, Polygon *pOut=0, bool instantReturnIfFind=true);
std::vector<CrossResult> findCrossPointsHorizontalOnPolygon(PolygonRef points, int y);
void findCrossPointsHorizontal(Polygon &out, Polygons &polygons, int y);
std::vector<CrossResult> findCrossPointsVerticalOnPolygon(PolygonRef points, int x);
void findCrossPointsVertical(Polygon &out, Polygons &polygons, int x);

bool isInsideHorizontal(PolygonRef polygon, Point p);
bool isPointInsideInPolygon(Polygons *pOutline, Point p);
Point findClosestPoint(PolygonRef polygon, Point p);
Point findClosestPoint(Polygons &polygons, Point p);
Point getInsidePoint(Polygon *pCrossPoints, Point p, Point prevP);

Point getAdvancePoint(Polygons outline, Point p, Point prevP);
void makePolygon(Polygon &out, int64_t x, std::vector<Point> &yarray);
void optimizePolyline(Polygons &optimizedPolygons, Polygons &polygons);
void generateZigzagPattern(Polygon &poly, Polygons &outline, int64_t x, int64_t y0, int64_t y1, int lineSpacing, int extrusionWidth, int verticalStride, PointMatrix &matrix, ZIGZAG_ALIGN align = ZIGZAGALIGN_CENTER);
void generateZigzagPattern1(Polygons &patterns, Polygons &outline, int64_t x, int64_t y0, int64_t y1, int lineSpacing, int extrusionWidth, int verticalStride, PointMatrix &matrix, ZIGZAG_ALIGN align);
void extendAABB(AABB &boundary, int offset);
void makeAABBOutline(PolygonRef outPoly, AABB boundary);
vector< vector<int64_t> > makeCutlist(AABB &boundary, Polygons &outline, int lineSpacing, int offset=0);
vector< vector<Point> > makeCutlistXY(AABB boundary, Polygons &outline, int lineSpacing, int offset=0);
void makeTriCutList(TriCutList &outList, AABB boundary, Polygons &outline, int lineSpacing, int offset=0);
void optimizePolylinePath(Polygons &optPolygons, Polygons &polygons, Point startPoint);
Polygons makeHeadCleaningPath(int minX, int maxX, int y, int count);
int readjustFanspeed(int fanSpeed, int fanThreshold);

Vector2 RotateVector(float radian, Vector2 p);
Vector2 FindPerpendicularVector(Point a, Point b, float rad, int distanceInMillimeter);
int FindIndexOfPointFromEnd(PolygonRef polygon, int distance);
int FindIndexOfPointFromBegin(PolygonRef polygon, int distance);
void TrimPolygonFromBegin(PolygonRef polygon, int distance);
void TrimPolygonFromEnd(PolygonRef polygon, int distance);
void ReOrderPolygon(PolygonRef polygon, int startIdx);
void MoveStartPosition(PolygonRef polygon, int offset);

void RearrangeStartPosition(Polygons &polygons, float angleInDeg);
void RearrangeOutlineStartPosition(SliceVolumeStorage& storage, float angleInDeg);

Polygons generateBedEdgeLine(int32_t distance, Point startP = Point(0, 0));
Polygons generateHomeLine(Polygons &polygons, int32_t maxDistance, Point from=Point(0, 0));
Polygons genInitPath(ConfigSettings& config, int32_t distance, Point startP = Point(0, 0), Polygons *pPolygons = NULL);

void writePolygon(const char *fileName, PolygonRef polygon,  const char *comment=NULL);
void writeOptions(GCodeExport &gcode, ConfigSettings &config, const char *fileName=NULL);

template< typename T, typename From >
const T StringCast(const From& from)
{
	std::stringstream strStream;
	strStream << from;

	T result;
	strStream >> result;

	if (strStream.fail() && !strStream.eof()) return 0;

	return result;
}
#endif//UTIL_H