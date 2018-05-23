/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef GCODEEXPORT_H
#define GCODEEXPORT_H

#include <stdio.h>

#include "settings.h"
#include "comb.h"
#include "utils/intpoint.h"
#include "utils/polygon.h"
#include "timeEstimate.h"
#include "GCodeControl.h"

#ifdef WIN32
#ifndef off64_t
#define off64_t __int64
#define lseek64 _lseeki64
#define fstat64 _fstati64
#define ftell64 _ftelli64
#define stat64 _stati64
#define ftello64 _ftelli64
#define fseeko64 _fseeki64
#endif

#endif

//The GCodeExport class writes the actual GCode. This is the only class that knows how GCode looks and feels.
//  Any customizations on GCodes flavors are done in this class.
class GCodeExport
{
private:
    FILE* f;
    double extrusionAmount;
    double extrusionPerMM;
    double retractionAmount;
    int retractionZHop;
    double extruderSwitchRetraction;
    double minimalExtrusionBeforeRetraction;
    double extrusionAmountAtPreviousRetraction;
    Point3 currentPosition;
    Point extruderOffset[MAX_EXTRUDERS];
    int currentSpeed, retractionSpeed;
    int zPos;
    bool isRetracted;
    int extruderNr;
    int currentFanSpeed;
    int flavor;
    
    double totalFilament[MAX_EXTRUDERS];
    double totalPrintTime;
    TimeEstimateCalculator estimateCalculator;

	int retractRestoreCorrection;
public:
    
    GCodeExport();
    
    ~GCodeExport();
    
    void replaceTagInStart(const char* tag, const char* replaceValue);
    
    void setExtruderOffset(int id, Point p);
    
    void setFlavor(int flavor);
    int getFlavor();
    
    void setFilename(const char* filename);
    
    bool isOpened();
    
    void setExtrusion(int layerThickness, int filamentDiameter, int flow);
    
	void setRetractionSettings(int retractionAmount, int retractionSpeed, int extruderSwitchRetraction, int minimalExtrusionBeforeRetraction, int zHop, int restoreCorrection);
    
    void setZ(int z);
    
    Point getPositionXY();
    
    int getPositionZ();

    int getExtruderNr();
    
    double getTotalFilamentUsed(int e);

    double getTotalPrintTime();
    void updateTotalPrintTime();
    
    void writeComment(const char* comment, ...);

    void writeLine(const char* line, ...);
    
    void resetExtrusionValue();
    
    void writeDelay(double timeAmount);
    
    void writeMove(Point p, int speed, int lineWidth);
    
    void writeRetraction();	

    void switchExtruder(int newExtruder);
    
    void writeCode(const char* str);
    
    void writeFanCommand(int speed);
    
	void finalize(int maxObjectHeight, int moveSpeed, const char* endCode, bool retraction = true, int extLift=1000);

    int getFileSize();
    void tellFileSize();

	//cstyle add new methods	
	void writeWipeExtrusion(Point p, double extrude);
	bool isNeedLayerRetraction(){ if (zPos != currentPosition.z) return true; return false;}
	void SetRetraction(bool retraction);
	int SetRetractionRestoreCorrection(int value){ retractRestoreCorrection = value; }
	int GetRectractionRestoreCorrection(){return retractRestoreCorrection;}
};

//The GCodePathConfig is the configuration for moves/extrusion actions. This defines at which width the line is printed and at which speed.
class GCodePathConfig
{
public:
    int speed;
    int lineWidth;
    const char* name;
    bool spiralize;
    
    GCodePathConfig() : speed(0), lineWidth(0), name(NULL), spiralize(false) {}
    GCodePathConfig(int speed, int lineWidth, const char* name) : speed(speed), lineWidth(lineWidth), name(name), spiralize(false) {}
    
    void setData(int speed, int lineWidth, const char* name)
    {
        this->speed = speed;
        this->lineWidth = lineWidth;
        this->name = name;
    }
};

class GCodePath
{
public:
    GCodePathConfig* config;
    bool retract;
    int extruder;
    vector<Point> points;
    bool done;//Path is finished, no more moves should be added, and a new path should be started instead of any appending done to this one.
	bool endRetract;
	//cstyle 
	double extrusionAmount;
	Point  wipePoint;

	//cstyle adjust extrusion amount for end of out wall
	int reduceStartIndex;//startIndex
	float reduceRate;//0~1
};

class GCodeControl;
//The GCodePlanner class stores multiple moves that are planned.
// It facilitates the combing to keep the head inside the print.
// It also keeps track of the print time estimate for this planning so speed adjustments can be made for the minimal-layer-time.
class GCodePlanner
{
private:
    GCodeExport& gcode;
    
    Point lastPosition;
    vector<GCodePath> paths;
    Comb* comb;
    
    GCodePathConfig travelConfig;
	GCodePathConfig travelConfig1;

    int extrudeSpeedFactor;
    int travelSpeedFactor;
    int currentExtruder;
    int retractionMinimalDistance;
    bool forceRetraction;
    bool alwaysRetract;
    double extraTime;
    double totalPrintTime;

	//cstyle
	int mStartSpeed;
private:
    
    void forceNewPathStart();
public:
	GCodePath* getLatestPathWithConfig(GCodePathConfig* config);
	GCodePlanner(GCodeExport& gcode, int travelSpeed, int retractionMinimalDistance, int startSpeed = 0);
    ~GCodePlanner();
    
    bool setExtruder(int extruder)
    {
        if (extruder == currentExtruder)
            return false;
        currentExtruder = extruder;
        return true;
    }
    
    int getExtruder()
    {
        return currentExtruder;
    }

    void setCombBoundary(Polygons* polygons)
    {
        if (comb)
            delete comb;
        if (polygons)
            comb = new Comb(*polygons);
        else
            comb = NULL;
    }
    
    void setAlwaysRetract(bool alwaysRetract)
    {
        this->alwaysRetract = alwaysRetract;
    }
    
    void forceRetract()
    {
        forceRetraction = true;
    }
    
    void setExtrudeSpeedFactor(int speedFactor)
    {
        if (speedFactor < 1) speedFactor = 1;
        this->extrudeSpeedFactor = speedFactor;
    }
    int getExtrudeSpeedFactor()
    {
        return this->extrudeSpeedFactor;
    }
    void setTravelSpeedFactor(int speedFactor)
    {
        if (speedFactor < 1) speedFactor = 1;
        this->travelSpeedFactor = speedFactor;
    }
    int getTravelSpeedFactor()
    {
        return this->travelSpeedFactor;
    }
    
	void addTravel(Point p, int speed=0);
	void forceAddTravel(Point p, int speed=0);
	void addExtrusionMove(Point p, GCodePathConfig* config, bool endRetraction = false);
	void SetEndReduce(GCodePathConfig* config, int endReduceIndex, float reduceRate = 1.0f);
    void moveInsideCombBoundary(int distance);

	void addPolygon(PolygonRef polygon, int startIdx, GCodePathConfig* config, GCodeControl ctrl=GCodeControl());// bool closeLine = true, bool forceRetract = false, bool endReduce = false, bool travelSpeedMatching = false, int endGap = 0);// , bool extrudeMoveAtBeginning = false);
	//void addPolygonsByOptimizer(Polygons& polygons, GCodePathConfig* config, bool keepPointOrder = false, Point *pStartPoint = NULL, bool endRetraction = false);// , bool extrudeMoveToStartPoint = false);
	//void addPolygonsByOptimizer(Polygons& polygons, GCodePathConfig* config, bool keepPointOrder=false, bool endRetraction=false);
	void addPolygonsByOptimizer(Polygons& polygons, GCodePathConfig* config, GCodeControl ctrl = GCodeControl());// , bool endRetraction = false, bool keepPointOrder = false, bool endReduce = false, bool travelSpeedMatching = false, int endGap = 0);
    void forceMinimalLayerTime(double minTime, int minimalSpeed);
    
    void writeGCode(bool liftHeadIfNeeded, int layerThickness, bool lockRetraction=false);

	void addWipeExtrusion(Point p, double extrusionAmount);
};

#endif//GCODEEXPORT_H
