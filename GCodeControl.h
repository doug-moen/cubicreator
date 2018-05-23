#ifndef GCODECONTROL_H
#define GCODECONTROL_H

class PolygonRef;
//cstyle
class GCodeControl
{
public:
	bool endRetraction;
	int retractionLimitDistance;
	bool keepPointOrder;

	bool travelSpeedMatching;
	int endGap;	

	//end reduce
	bool endReduce;
	int reduceDist;	
	float reduceRate;

	int outwallStartOffset;

	int startEntryAngle;
	int startEntryTravelDistance;

	bool useSameStartPosition;

	GCodeControl() :
		endRetraction(false),
		retractionLimitDistance(30000),//limit retraction and endreduce when polygon length is smaller than 30mm
		keepPointOrder(false),
		travelSpeedMatching(false),
		endGap(0),
		endReduce(false),
		reduceDist(0),		
		reduceRate(1.0f),
		outwallStartOffset(0),
		startEntryAngle(0),
		startEntryTravelDistance(0),
		useSameStartPosition(false)
	{
	}

	void ApplyEndGap(PolygonRef polygon);
	int ApplyEndReduce(PolygonRef polygon);

};

#endif//GCODECONTROL_H