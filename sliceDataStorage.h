/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#ifndef SLICE_DATA_STORAGE_H
#define SLICE_DATA_STORAGE_H

#include "utils/intpoint.h"
#include "utils/polygon.h"

/*
SliceData
+ Layers[]
  + LayerParts[]
    + OutlinePolygons[]
    + Insets[]
      + Polygons[]
    + SkinPolygons[]
*/

class SliceLayerPart
{
public:
    AABB boundaryBox;
    Polygons outline;
    Polygons combBoundery;
    vector<Polygons> insets;
    Polygons skinOutline;
    Polygons sparseOutline;
    int bridgeAngle;
};

class SliceLayer
{
public:
    int z;
    vector<SliceLayerPart> parts;
};

/******************/
class SupportPoint
{
public:
    int32_t z;
    double cosAngle;
    
    SupportPoint(int32_t z, double cosAngle) : z(z), cosAngle(cosAngle) {}
};
class SupportStorage
{
public:
    bool generated;
    int angle;
    bool everywhere;
    int XYDistance;
    int ZDistance;
    
    Point gridOffset;
    int32_t gridScale;
    int32_t gridWidth, gridHeight;
    vector<SupportPoint>* grid;
   	SupportStorage(){grid = NULL;}
	  ~SupportStorage(){if(grid) delete [] grid;}
};
/******************/

class SliceVolumeStorage
{
public:
    vector<SliceLayer> layers;
	
	void RemoveEmptyLayer()
	{
		for (int i = 0; i < layers.size(); i++)
		{
			if (layers[i].parts.size() == 0)
			{
				layers.resize(i);
				//layers.shrink_to_fit();
			}
		}
	}	

	Polygons GetPolygons()
	{
		Polygons polys;
		for (int i = 0; i < layers.size(); i++)
		{
			for (int k = 0; k < layers[i].parts.size(); k++)
			{
				polys.add(layers[i].parts[k].insets[0]);
				//layers.shrink_to_fit();
			}
		}
		return polys;
	}
};

class SliceDataStorage
{
public:
    Point3 modelSize, modelMin, modelMax;
    Polygons skirt;
    Polygons raftOutline;               //Storage for the outline of the raft. Will be filled with lines when the GCode is generated.
	Polygons raftLiftOutline;//cstyle
    vector<Polygons> oozeShield;        //oozeShield per layer
    vector<SliceVolumeStorage> volumes;
    
    SupportStorage support;
    Polygons wipeTower;
    Point wipePoint;
	//cstyle 
	Point wipePoint1;
	Polygons boundaryBox;
};

#endif//SLICE_DATA_STORAGE_H
