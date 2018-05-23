/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <stdarg.h>

#include "gcodeExport.h"
#include "pathOrderOptimizer.h"
#include "timeEstimate.h"
#include "settings.h"
#include "utils/logoutput.h"
#include "utils/util.h"

#if defined(__APPLE__) && defined(__MACH__)
//On MacOS the file offset functions are always 64bit.
#define off64_t off_t
#define ftello64 ftello
#define fseeko64 fseeko
#else
#define off64_t int
#define ftello64 ftell
#define fseeko64 fseek
#endif

GCodeExport::GCodeExport()
: currentPosition(0,0,0),
retractRestoreCorrection(0)
{
    extrusionAmount = 0;
    extrusionPerMM = 0;
    retractionAmount = 4.5;
    minimalExtrusionBeforeRetraction = 0.0;
    extrusionAmountAtPreviousRetraction = -10000;
    extruderSwitchRetraction = 14.5;
    extruderNr = 0;
    currentFanSpeed = -1;
    
    totalPrintTime = 0.0;
    for(unsigned int e=0; e<MAX_EXTRUDERS; e++)
        totalFilament[e] = 0.0;
    
    currentSpeed = 0;
    retractionSpeed = 45;
    isRetracted = true;
    memset(extruderOffset, 0, sizeof(extruderOffset));
    f = stdout;	
}

GCodeExport::~GCodeExport()
{
    if (f && f != stdout)
        fclose(f);
}

void GCodeExport::replaceTagInStart(const char* tag, const char* replaceValue)
{
    if (f == stdout)
    {
        log("Replace:%s:%s\n", tag, replaceValue);
        return;
    }
    off64_t oldPos = ftello64(f);
    
    char buffer[1024];
    fseeko64(f, 0, SEEK_SET);
    fread(buffer, 1024, 1, f);
    
    char* c = strstr(buffer, tag);
	if (c)
	{
		memset(c, ' ', strlen(tag));
		if (c) memcpy(c, replaceValue, strlen(replaceValue));
	}
    
    fseeko64(f, 0, SEEK_SET);
    fwrite(buffer, 1024, 1, f);
    
    fseeko64(f, oldPos, SEEK_SET);
}

void GCodeExport::setExtruderOffset(int id, Point p)
{
    extruderOffset[id] = p;
}

void GCodeExport::setFlavor(int flavor)
{
    this->flavor = flavor;
}
int GCodeExport::getFlavor()
{
    return this->flavor;
}

void GCodeExport::setFilename(const char* filename)
{
    f = fopen(filename, "w+");
}

bool GCodeExport::isOpened()
{
    return f != NULL;
}

void GCodeExport::setExtrusion(int layerThickness, int filamentDiameter, int flow)
{
    double filamentArea = M_PI * (double(filamentDiameter) / 1000.0 / 2.0) * (double(filamentDiameter) / 1000.0 / 2.0);
    if (flavor == GCODE_FLAVOR_ULTIGCODE)//UltiGCode uses volume extrusion as E value, and thus does not need the filamentArea in the mix.
        extrusionPerMM = double(layerThickness) / 1000.0;
    else
        extrusionPerMM = double(layerThickness) / 1000.0 / filamentArea * double(flow) / 100.0;
}

void GCodeExport::setRetractionSettings(int retractionAmount, int retractionSpeed, int extruderSwitchRetraction, int minimalExtrusionBeforeRetraction, int zHop, int restoreCorrection)
{
    this->retractionAmount = double(retractionAmount) / 1000.0;
    this->retractionSpeed = retractionSpeed;
    this->extruderSwitchRetraction = double(extruderSwitchRetraction) / 1000.0;
    this->minimalExtrusionBeforeRetraction = double(minimalExtrusionBeforeRetraction) / 1000.0;
    this->retractionZHop = zHop;
	retractRestoreCorrection = restoreCorrection;
}

void GCodeExport::setZ(int z)
{
    this->zPos = z;
}

Point GCodeExport::getPositionXY()
{
    return Point(currentPosition.x, currentPosition.y);
}

int GCodeExport::getPositionZ()
{
    return currentPosition.z;
}

int GCodeExport::getExtruderNr()
{
    return extruderNr;
}

double GCodeExport::getTotalFilamentUsed(int e)
{
    if (e == extruderNr)
        return totalFilament[e] + extrusionAmount;
    return totalFilament[e];
}

double GCodeExport::getTotalPrintTime()
{
    return totalPrintTime;
}

void GCodeExport::updateTotalPrintTime()
{
    totalPrintTime += estimateCalculator.calculate();
    estimateCalculator.reset();
}

void GCodeExport::writeComment(const char* comment, ...)
{
    va_list args;
    va_start(args, comment);
    fprintf(f, ";");
    vfprintf(f, comment, args);
    fprintf(f, "\n");
    va_end(args);
}

void GCodeExport::writeLine(const char* line, ...)
{
    va_list args;
    va_start(args, line);
    vfprintf(f, line, args);
    fprintf(f, "\n");
    va_end(args);
}

void GCodeExport::resetExtrusionValue()
{
    if (extrusionAmount != 0.0)
    {
        fprintf(f, "G92 E0\n");
        totalFilament[extruderNr] += extrusionAmount;
        extrusionAmountAtPreviousRetraction -= extrusionAmount;
        extrusionAmount = 0.0;
    }
}

void GCodeExport::writeDelay(double timeAmount)
{
    fprintf(f, "G4 P%d\n", int(timeAmount * 1000));
    totalPrintTime += timeAmount;
}

//org
void GCodeExport::writeMove(Point p, int speed, int lineWidth)
{
	if (lineWidth != 0)
	{
		Point diff = p - getPositionXY();
		if (isRetracted)
		{
			if (retractionZHop > 0)
				fprintf(f, "G1 Z%0.3f\n", float(currentPosition.z) / 1000);
			if (flavor == GCODE_FLAVOR_ULTIGCODE)
			{
				fprintf(f, "G11\n");
			}
			else
			{
				fprintf(f, "G1 F%i E%0.5lf\n", retractionSpeed * 60, extrusionAmount);
				currentSpeed = retractionSpeed;
				estimateCalculator.plan(TimeEstimateCalculator::Position(double(p.X) / 1000.0, (p.Y) / 1000.0, double(zPos) / 1000.0, extrusionAmount), currentSpeed);
			}
			if (extrusionAmount > 10000.0) //According to https://github.com/Ultimaker/CuraEngine/issues/14 having more then 21m of extrusion causes inaccuracies. So reset it every 10m, just to be sure.				
				resetExtrusionValue();
			isRetracted = false;
		}
		extrusionAmount += extrusionPerMM * double(lineWidth) / 1000.0 * vSizeMM(diff);
		fprintf(f, "G1");
	}
	else
	{
		//cstyle z single move test
		//if (zPos != currentPosition.z)
		//{
		//	fprintf(f, "G0");
		//	if (currentSpeed != speed)
		//	{
		//		fprintf(f, " F%i", speed * 60);
		//		currentSpeed = speed;
		//	}
		//	fprintf(f, " Z%0.3f\n", float(zPos) / 1000);
		//}
		//cstyle End z single move test

		fprintf(f, "G0");
	}
	
	if (currentSpeed != speed)
	{
		fprintf(f, " F%i", speed * 60);
		currentSpeed = speed;
	}

	fprintf(f, " X%0.3f Y%0.3f", float(p.X - extruderOffset[extruderNr].X) / 1000, float(p.Y - extruderOffset[extruderNr].Y) / 1000);
	if (zPos != currentPosition.z)//org
	//if (lineWidth != 0 && zPos != currentPosition.z)//cstyle z single move test
	    fprintf(f, " Z%0.3f", float(zPos)/1000);
	if (lineWidth != 0)
		fprintf(f, " E%0.5lf", extrusionAmount);
	fprintf(f, "\n");

	currentPosition = Point3(p.X, p.Y, zPos);
	estimateCalculator.plan(TimeEstimateCalculator::Position(double(currentPosition.x) / 1000.0, (currentPosition.y) / 1000.0, double(currentPosition.z) / 1000.0, extrusionAmount), currentSpeed);
}

void GCodeExport::writeRetraction()
{
    if (retractionAmount > 0 && !isRetracted && extrusionAmountAtPreviousRetraction + minimalExtrusionBeforeRetraction < extrusionAmount)
    {
        if (flavor == GCODE_FLAVOR_ULTIGCODE)
        {
            fprintf(f, "G10\n");
        }
		else
		{
			//cstyle  - correction value
			//double corrVal = (double)retractRestoreCorrection/1000.0f;
			//extrusionAmount += corrVal;
            //fprintf(f, "G1 F%i E%0.5lf\n", retractionSpeed * 60, extrusionAmount - retractionAmount+corrVal);
			fprintf(f, "G1 F%i E%0.5lf\n", retractionSpeed * 60, extrusionAmount - retractionAmount);
            currentSpeed = retractionSpeed;
            estimateCalculator.plan(TimeEstimateCalculator::Position(double(currentPosition.x) / 1000.0, (currentPosition.y) / 1000.0, double(currentPosition.z) / 1000.0, extrusionAmount - retractionAmount), currentSpeed);
        }
        if (retractionZHop > 0)
            fprintf(f, "G1 Z%0.2f\n", float(currentPosition.z + retractionZHop)/1000);
        extrusionAmountAtPreviousRetraction = extrusionAmount;
        isRetracted = true;
    }
}

void GCodeExport::switchExtruder(int newExtruder)
{
    if (extruderNr == newExtruder)
        return;
    
    resetExtrusionValue();
    extruderNr = newExtruder;

    if (flavor == GCODE_FLAVOR_ULTIGCODE)
    {
        fprintf(f, "G10 S1\n");
    }else{
        fprintf(f, "G1 F%i E%0.4lf\n", retractionSpeed * 60, extrusionAmount - extruderSwitchRetraction);
        currentSpeed = retractionSpeed;
    }
    isRetracted = true;
    if (flavor == GCODE_FLAVOR_MAKERBOT)
        fprintf(f, "M135 T%i\n", extruderNr);
    else
        fprintf(f, "T%i\n", extruderNr);
}

void GCodeExport::writeCode(const char* str)
{
    fprintf(f, "%s\n", str);
}

//speed range:0~100
void GCodeExport::writeFanCommand(int speed)
{
    if (currentFanSpeed == speed)
        return;
    if (speed > 0)
    {
        if (flavor == GCODE_FLAVOR_MAKERBOT)
            fprintf(f, "M126 T0 ; value = %d\n", speed * 255 / 100);
        else
            fprintf(f, "M106 S%d\n", speed * 255 / 100);
    }
    else
    {
        if (flavor == GCODE_FLAVOR_MAKERBOT)
            fprintf(f, "M127 T0\n");
        else
            fprintf(f, "M107\n");
    }
    currentFanSpeed = speed;
}

int GCodeExport::getFileSize(){
    return ftell(f);
}
void GCodeExport::tellFileSize() {
    float fsize = (float) ftell(f);
    if(fsize > 1024*1024) {
        fsize /= 1024.0*1024.0;
        log("Wrote %5.1f MB.\n",fsize);
    }
    if(fsize > 1024) {
        fsize /= 1024.0;
        log("Wrote %5.1f kilobytes.\n",fsize);
    }
}

void GCodeExport::finalize(int maxObjectHeight, int moveSpeed, const char* endCode, bool retraction, int extLift)
{
    writeFanCommand(0);    
	if (retraction)
		writeRetraction();
	setZ(maxObjectHeight + extLift);
    writeMove(getPositionXY(), moveSpeed, 0);
    writeCode(endCode);
    log("Print time: %d\n", int(getTotalPrintTime()));
    log("Filament: %d\n", int(getTotalFilamentUsed(0)));
    log("Filament2: %d\n", int(getTotalFilamentUsed(1)));
    
    if (getFlavor() == GCODE_FLAVOR_ULTIGCODE)
    {
        char numberString[16];
        sprintf(numberString, "%d", int(getTotalPrintTime()));
        replaceTagInStart("<__TIME__>", numberString);
        sprintf(numberString, "%d", int(getTotalFilamentUsed(0)));
        replaceTagInStart("<FILAMENT>", numberString);
        sprintf(numberString, "%d", int(getTotalFilamentUsed(1)));
        replaceTagInStart("<FILAMEN2>", numberString);
    } 
	else if (getFlavor() == GCODE_FLAVOR_CUBICREATOR)
	{
		char numberString[16];
		sprintf(numberString, "%d", int(getTotalPrintTime()));
		replaceTagInStart("<__TIME__>", numberString);
		memset(numberString, 0x00, 16);
		sprintf(numberString, "%d", int(getFileSize()));
		replaceTagInStart("<__SIZE__>", numberString);
		sprintf(numberString, "%d", int(getTotalFilamentUsed(0)));
		replaceTagInStart("<FILAMENT>", numberString);
		sprintf(numberString, "%d", int(getTotalFilamentUsed(1)));
		replaceTagInStart("<FILAMENT2>", numberString);
	}

}

void GCodeExport::writeWipeExtrusion(Point p, double amount)
{
	int		extrusionSpeed	= 200;
	double	retractionAmount= -16.5f;
	int		travelSpeed		= 6000;

	if (flavor == GCODE_FLAVOR_ULTIGCODE)
    {
        fprintf(f, "G10\n");
    } else
	{
		float x = float(p.X - extruderOffset[extruderNr].X)/1000;
		float y = float(p.Y - extruderOffset[extruderNr].Y)/1000;

		fprintf(f, "G91\n");//relative position
		//fprintf(f, "G1 F300\n");//fprintf(f, "G1 E-1 F300\n");
		fprintf(f, "G1 Z+0.5 F600\n");//fprintf(f, "G1 Z+0.5 E-5 F600\n");
		fprintf(f, "G90\n");//absolute position

		//fprintf(f, "G0 F%i X%0.2f Y%0.2f\n", travelSpeed, x, y);//travels to wipe point.
		fprintf(f, "G0 F%i X%0.2f\n", travelSpeed, x);//travels to wipe point.

		extrusionAmount += (amount/1000);
        //fprintf(f, "G1 F%i X%0.2f Y%0.2f E%0.5lf\n", 200, x, y, extrusionAmount);//retractionSpeed * 60, extrusionAmount);
		fprintf(f, "G1 F%i E%0.5lf\n", extrusionSpeed, extrusionAmount);
		resetExtrusionValue();
		//fprintf(f, "G92 E-5.0\n");//retractionSpeed * 60, extrusionAmount);
		//fprintf(f, "G0 F%i X%0.2f Y%0.2f\n", travelSpeed, x, y);
		//fprintf(f, "G1 F%i E%0.5lf\n", extrusionSpeed, retractionAmount);

		fprintf(f, "G91\n");//relative position		
		fprintf(f, "G1 Z-0.5 F600\n");//fprintf(f, "G1 Z+0.5 E-5 F600\n");
		fprintf(f, "G90\n");//absolute position

        currentSpeed = retractionSpeed;
        estimateCalculator.plan(TimeEstimateCalculator::Position(double(currentPosition.x) / 1000.0, (currentPosition.y) / 1000.0, double(currentPosition.z) / 1000.0, extrusionAmount - retractionAmount), currentSpeed);
    }
}

void GCodeExport::SetRetraction(bool retraction)
{
	isRetracted = retraction;
}

GCodePath* GCodePlanner::getLatestPathWithConfig(GCodePathConfig* config)
{
    if (paths.size() > 0 && paths[paths.size()-1].config == config && !paths[paths.size()-1].done)
        return &paths[paths.size()-1];
    paths.push_back(GCodePath());
    GCodePath* ret = &paths[paths.size()-1];
    ret->retract = false;
    ret->config = config;
    ret->extruder = currentExtruder;
    ret->done = false;
	
	//cstyle test code for cleaning header.
	ret->extrusionAmount = 0;
	ret->wipePoint = Point(0, 0);
	ret->endRetract = false;
	//cstyle
	ret->reduceStartIndex = -1;
	ret->reduceRate = 1.0f;

    return ret;
}
void GCodePlanner::forceNewPathStart()
{
    if (paths.size() > 0)
        paths[paths.size()-1].done = true;
}

GCodePlanner::GCodePlanner(GCodeExport& gcode, int travelSpeed, int retractionMinimalDistance, int startSpeed): 
gcode(gcode), 
travelConfig(travelSpeed, 0, "travel"), 
travelConfig1(20, 0, "travel1"),
mStartSpeed(startSpeed)
{
    lastPosition = gcode.getPositionXY();
    comb = NULL;
    extrudeSpeedFactor = 100;
	travelSpeedFactor = 100;
    extraTime = 0.0;
    totalPrintTime = 0.0;
    forceRetraction = false;
    alwaysRetract = false;
    currentExtruder = gcode.getExtruderNr();
    this->retractionMinimalDistance = retractionMinimalDistance;
}
GCodePlanner::~GCodePlanner()
{
    if (comb)
        delete comb;
}

void GCodePlanner::addTravel(Point p, int speed)
{
	//cstyle
	if (speed != 0)
		travelConfig1.speed = speed;
	//GCodePath* path = getLatestPathWithConfig(&travelConfig); org
    GCodePath* path = getLatestPathWithConfig(speed!=0?&travelConfig1:&travelConfig);
    if (forceRetraction)
    {
        if (!shorterThen(lastPosition - p, retractionMinimalDistance))
        {
            path->retract = true;
        }
        forceRetraction = false;
    }else if (comb != NULL)
    {
        vector<Point> pointList;
        if (comb->calc(lastPosition, p, pointList))
        {
            for(unsigned int n=0; n<pointList.size(); n++)
            {
                path->points.push_back(pointList[n]);
            }
        }else{
            if (!shorterThen(lastPosition - p, retractionMinimalDistance))
                path->retract = true;
        }
    }else if (alwaysRetract)
    {
        if (!shorterThen(lastPosition - p, retractionMinimalDistance))
            path->retract = true;
    }
    path->points.push_back(p);
    lastPosition = p;
}

//cstyle custom function
void GCodePlanner::forceAddTravel(Point p, int speed)
{
	if (speed!=0)
		travelConfig1.speed = speed;	
	GCodePath* path = getLatestPathWithConfig(speed != 0 ? &travelConfig1 : &travelConfig);
	path->points.push_back(p);
	lastPosition = p;
}

void GCodePlanner::addExtrusionMove(Point p, GCodePathConfig* config, bool endRetraction)
{
    //getLatestPathWithConfig(config)->points.push_back(p);//org
	GCodePath *pPath = getLatestPathWithConfig(config);
	pPath->points.push_back(p);
	pPath->endRetract = endRetraction;
    lastPosition = p;
}

void GCodePlanner::SetEndReduce(GCodePathConfig* config, int endReduceIndex, float reduceRate)
{
	GCodePath *pPath = getLatestPathWithConfig(config);
	if (endReduceIndex!=-1)
	{
		pPath->reduceStartIndex = endReduceIndex;
		pPath->reduceRate = reduceRate;
	}
}

void GCodePlanner::moveInsideCombBoundary(int distance)
{
    if (!comb || comb->checkInside(lastPosition)) return;
    Point p = lastPosition;
    if (comb->moveInside(&p, distance))
    {
        //Move inside again, so we move out of tight 90deg corners
        comb->moveInside(&p, distance);
        if (comb->checkInside(p))
        {
            addTravel(p);
            //Make sure the that any retraction happens after this move, not before it by starting a new move path.
            forceNewPathStart();
        }
    }
}

void GCodePlanner::addPolygon(PolygonRef polygon, int startIdx, GCodePathConfig* config, GCodeControl ctrl)
{
	if (polygon.size() == 0)
		return;

	if (startIdx != 0)
		ReOrderPolygon(polygon, startIdx);	

	//move start position. 
	if (ctrl.outwallStartOffset != 0)	
		MoveStartPosition(polygon, ctrl.outwallStartOffset);
	
    Point p0 = polygon[0];

	//addTravel(p0);org
	if (ctrl.startEntryTravelDistance > 0 && polygon.size()>2)
	{
		Vector2 v2 = FindPerpendicularVector(polygon[0], polygon[1], DEGTORAD((float)ctrl.startEntryAngle), ctrl.startEntryTravelDistance);
		Point perpendicularP = Point((int)v2.X, (int)v2.Y);//FindPerpendicularVector(polygon[0], polygon[1], DEGTORAD((float)ctrl.startEntryAngle), ctrl.startEntryTravelDistance);
		forceAddTravel(perpendicularP, config->speed);
	}
	addTravel(p0, ctrl.travelSpeedMatching?config->speed:0);	
	
	//end gap
	ctrl.ApplyEndGap(polygon);	
    for(unsigned int i=1; i<polygon.size(); i++)
    {
		int index = (0 + i) % polygon.size();
        Point p1 = polygon[index];

		if (ctrl.endGap>0 && i == polygon.size() - 1)
			addExtrusionMove(p1, config, ctrl.endRetraction);
		else
			addExtrusionMove(p1, config, false);		
    }

	//this makes loop line. connect the start point with the end point.
	if (ctrl.endGap == 0 && (!ctrl.keepPointOrder || ctrl.useSameStartPosition) && polygon.size() > 2)
	{
		Point last = polygon[polygon.size() - 1];
		if (last.X!=polygon[0].X || last.Y!=polygon[0].Y)
			addExtrusionMove(polygon[0], config, ctrl.endRetraction);//cstyle org
	}

	//cstyle this endReduce is for only outwall
	int reduceIndex = ctrl.ApplyEndReduce(polygon);
	if (ctrl.endReduce && reduceIndex>0)
		SetEndReduce(config, reduceIndex, ctrl.reduceRate);

	if (ctrl.endGap > 0)
		lastPosition = polygon[0];//cstyle org		
}

void GCodePlanner::addPolygonsByOptimizer(Polygons& polygons, GCodePathConfig* config, GCodeControl ctrl)
{		
    PathOrderOptimizer orderOptimizer(lastPosition);
    for(unsigned int i=0;i<polygons.size();i++)
        orderOptimizer.addPolygon(polygons[i]);
	orderOptimizer.optimize(ctrl.keepPointOrder || ctrl.useSameStartPosition);
	
    for(unsigned int i=0;i<orderOptimizer.polyOrder.size();i++)
    {
        int nr = orderOptimizer.polyOrder[i];
		//cstyle endretraction		
		if (polygons[nr].polygonLength() < ctrl.retractionLimitDistance)
		{
			ctrl.endReduce = false;
			ctrl.endRetraction = false;
		}
		else
			ctrl.endRetraction = (ctrl.endRetraction && i == orderOptimizer.polyOrder.size() - 1) ? true : false;

		addPolygon(polygons[nr], orderOptimizer.polyStart[nr], config, ctrl);
    }	
}

void GCodePlanner::forceMinimalLayerTime(double minTime, int minimalSpeed)
{
    Point p0 = gcode.getPositionXY();
    double travelTime = 0.0;
    double extrudeTime = 0.0;
    for(unsigned int n=0; n<paths.size(); n++)
    {
        GCodePath* path = &paths[n];
        for(unsigned int i=0; i<path->points.size(); i++)
        {
            double thisTime = vSizeMM(p0 - path->points[i]) / double(path->config->speed);
            if (path->config->lineWidth != 0)
                extrudeTime += thisTime;
            else
                travelTime += thisTime;
            p0 = path->points[i];
        }
    }
    double totalTime = extrudeTime + travelTime;
    if (totalTime < minTime && extrudeTime > 0.0)
    {
        double minExtrudeTime = minTime - travelTime;
        if (minExtrudeTime < 1)
            minExtrudeTime = 1;
        double factor = extrudeTime / minExtrudeTime;
        for(unsigned int n=0; n<paths.size(); n++)
        {
            GCodePath* path = &paths[n];
            if (path->config->lineWidth == 0)
                continue;
            int speed = path->config->speed * factor;
            if (speed < minimalSpeed)
                factor = double(minimalSpeed) / double(path->config->speed);
        }
        
        //Only slow down with the minimal time if that will be slower then a factor already set. First layer slowdown also sets the speed factor.
        if (factor * 100 < getExtrudeSpeedFactor())
            setExtrudeSpeedFactor(factor * 100);
        else
            factor = getExtrudeSpeedFactor() / 100.0;
        
        if (minTime - (extrudeTime / factor) - travelTime > 0.1)
        {
            //TODO: Use up this extra time (circle around the print?)
            this->extraTime = minTime - (extrudeTime / factor) - travelTime;
        }
        this->totalPrintTime = (extrudeTime / factor) + travelTime;
    }else{
        this->totalPrintTime = totalTime;
    }
}

void GCodePlanner::writeGCode(bool liftHeadIfNeeded, int layerThickness, bool lockRetraction)
{
    GCodePathConfig* lastConfig = NULL;
    int extruder = gcode.getExtruderNr();

	gcode.SetRetraction(!lockRetraction);//cstyle test I don't know it work properly

    for(unsigned int n=0; n<paths.size(); n++)
    {
        GCodePath* path = &paths[n];
        if (extruder != path->extruder)
        {
            extruder = path->extruder;
            gcode.switchExtruder(extruder);

			//cstyle 
			if (path->extrusionAmount>0)
			{
				gcode.writeWipeExtrusion(path->wipePoint, path->extrusionAmount);
				path->extrusionAmount = 0;
			}

        }else if (path->retract && !lockRetraction)//cstyle just put a locking code
        {
            gcode.writeRetraction();
        }
        if (path->config != &travelConfig && lastConfig != path->config)
        {
            gcode.writeComment("TYPE:%s", path->config->name);
            lastConfig = path->config;
			//gcode.resetExtrusionValue();//cstyle test
        }

		//cstyle this is for initial printer work.
		//this process is needed only one time.
        int speed = path->config->speed;
		if (mStartSpeed > 0)
		{
			speed = mStartSpeed;
			mStartSpeed = 0;
		}
        
        if (path->config->lineWidth != 0)// Only apply the extrudeSpeedFactor to extrusion moves
            speed = speed * extrudeSpeedFactor / 100;
        else
            speed = speed * travelSpeedFactor / 100;
        
        if (path->points.size() == 1 && path->config != &travelConfig && shorterThen(gcode.getPositionXY() - path->points[0], path->config->lineWidth * 2))
        {
            //Check for lots of small moves and combine them into one large line
            Point p0 = path->points[0];
            unsigned int i = n + 1;
            while(i < paths.size() && paths[i].points.size() == 1 && shorterThen(p0 - paths[i].points[0], path->config->lineWidth * 2))
            {
                p0 = paths[i].points[0];
                i ++;
            }
            if (paths[i-1].config == &travelConfig)
                i --;
            if (i > n + 2)
            {
                p0 = gcode.getPositionXY();
                for(unsigned int x=n; x<i-1; x+=2)
                {
                    int64_t oldLen = vSize(p0 - paths[x].points[0]);
                    Point newPoint = (paths[x].points[0] + paths[x+1].points[0]) / 2;
                    int64_t newLen = vSize(gcode.getPositionXY() - newPoint);
                    if (newLen > 0)
                        gcode.writeMove(newPoint, speed, path->config->lineWidth * oldLen / newLen);
                    
                    p0 = paths[x+1].points[0];
                }
                gcode.writeMove(paths[i-1].points[0], speed, path->config->lineWidth);
                n = i - 1;
                continue;
            }
        }
        
        bool spiralize = path->config->spiralize;
        if (spiralize)
        {
            //Check if we are the last spiralize path in the list, if not, do not spiralize.
            for(unsigned int m=n+1; m<paths.size(); m++)
            {
                if (paths[m].config->spiralize)
                    spiralize = false;
            }
        }
        if (spiralize)
        {
            //If we need to spiralize then raise the head slowly by 1 layer as this path progresses.
            float totalLength = 0.0;
            int z = gcode.getPositionZ();
            Point p0 = gcode.getPositionXY();
            for(unsigned int i=0; i<path->points.size(); i++)
            {
                Point p1 = path->points[i];
                totalLength += vSizeMM(p0 - p1);
                p0 = p1;
            }
            
            float length = 0.0;
            p0 = gcode.getPositionXY();
            for(unsigned int i=0; i<path->points.size(); i++)
            {
                Point p1 = path->points[i];
                length += vSizeMM(p0 - p1);
                p0 = p1;
                gcode.setZ(z + layerThickness * length / totalLength);
                gcode.writeMove(path->points[i], speed, path->config->lineWidth);
            }
        }else{
            for(unsigned int i=0; i<path->points.size(); i++)
            {
                //gcode.writeMove(path->points[i], speed, path->config->lineWidth);//cstyle org
				if (path->reduceStartIndex!=-1 && i>=path->reduceStartIndex)
					//gcode.writeMove(path->points[i], speed, path->config->lineWidth*path->reduceRate);
					gcode.writeMove(path->points[i], speed*path->reduceRate, path->config->lineWidth);
				else
					gcode.writeMove(path->points[i], speed, path->config->lineWidth);
            }
        }

		if (path->endRetract && !lockRetraction)//cstyle just put a lock code
		{
			gcode.writeRetraction();
		}
    }
    
    gcode.updateTotalPrintTime();
    if (liftHeadIfNeeded && extraTime > 0.0)
    {
        gcode.writeComment("Small layer, adding delay of %f", extraTime);
        gcode.writeRetraction();
        gcode.setZ(gcode.getPositionZ() + 3000);
        gcode.writeMove(gcode.getPositionXY(), travelConfig.speed, 0);
        gcode.writeMove(gcode.getPositionXY() - Point(-20000, 0), travelConfig.speed, 0);
        gcode.writeDelay(extraTime);
    }
}

//cstyle
void GCodePlanner::addWipeExtrusion(Point p, double extrusionAmount)
{
	GCodePath* path = getLatestPathWithConfig(&travelConfig);
	if (path)
	{
		path->extrusionAmount = extrusionAmount;
		path->wipePoint = p;
		//path->points.push_back(p);
		lastPosition = p;
	}
}