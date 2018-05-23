/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include "raft.h"
#include "support.h"
#include "utils/util.h"

void generateRaft(SliceDataStorage& storage, int distance)
{
    for(unsigned int volumeIdx = 0; volumeIdx < storage.volumes.size(); volumeIdx++)
    {
        if (storage.volumes[volumeIdx].layers.size() < 1) continue;
        SliceLayer* layer = &storage.volumes[volumeIdx].layers[0];
        for(unsigned int i=0; i<layer->parts.size(); i++)
        {
            storage.raftOutline = storage.raftOutline.unionPolygons(layer->parts[i].outline.offset(distance));
        }
    }

    SupportPolyGenerator supportGenerator(storage.support, 0);
    storage.raftOutline = storage.raftOutline.unionPolygons(supportGenerator.polygons);
}

void generateRaft(SliceDataStorage& storage, int distance, int dashStride, int dashSize)
{
    for(unsigned int volumeIdx = 0; volumeIdx < storage.volumes.size(); volumeIdx++)
    {
        if (storage.volumes[volumeIdx].layers.size() < 1) continue;
        SliceLayer* layer = &storage.volumes[volumeIdx].layers[0];
        for(unsigned int i=0; i<layer->parts.size(); i++)
        {
			storage.raftOutline = storage.raftOutline.unionPolygons(layer->parts[i].outline);// .offset(distance));
			storage.raftLiftOutline = storage.raftLiftOutline.unionPolygons(layer->parts[i].outline);
        }		
    }
    SupportPolyGenerator supportGenerator(storage.support, 0);
	//cstyle
    storage.raftOutline = storage.raftOutline.unionPolygons(supportGenerator.polygons);	
	storage.raftOutline = storage.raftOutline.offset(distance);
}
