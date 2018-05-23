#ifndef FFF_PROCESSOR_H
#define FFF_PROCESSOR_H

//#include "utils/socket.h"
#include "utils/util.h"

#define GUI_CMD_REQUEST_MESH 0x01
#define GUI_CMD_SEND_POLYGONS 0x02

//FusedFilamentFabrication processor.
class fffProcessor
{
private:
    int maxObjectHeight;
    int fileNr;
    GCodeExport gcode;
    ConfigSettings& config;
    TimeKeeper timeKeeper;
    //ClientSocket guiSocket;

    GCodePathConfig skirtConfig;
    GCodePathConfig inset0Config;
    GCodePathConfig insetXConfig;
    GCodePathConfig fillConfig;
    GCodePathConfig supportConfig;
	GCodePathConfig homelineConfig;

	//cstyle
	bool drawHomeLine;
	std::string inputFileName;
public:
    fffProcessor(ConfigSettings& config)
		: config(config), drawHomeLine(true), inputFileName()
    {
        fileNr = 1;
        maxObjectHeight = 0;
    }
    
    void guiConnect(int portNr)
    {
        //guiSocket.connectTo("127.0.0.1", portNr);
    }
    
    void sendPolygonsToGui(const char* name, int layerNr, int32_t z, Polygons& polygons)
    {
        /*guiSocket.sendNr(GUI_CMD_SEND_POLYGONS);
        guiSocket.sendNr(polygons.size());
        guiSocket.sendNr(layerNr);
        guiSocket.sendNr(z);
        guiSocket.sendNr(strlen(name));
        guiSocket.sendAll(name, strlen(name));
        for(unsigned int n=0; n<polygons.size(); n++)
        {
            PolygonRef polygon = polygons[n];
            guiSocket.sendNr(polygon.size());
            guiSocket.sendAll(polygon.data(), polygon.size() * sizeof(Point));
        }*/
    }
    
    bool setTargetFile(const char* filename)
    {
        gcode.setFilename(filename);
        if (gcode.isOpened())
		{
            //gcode.writeComment("Generated with Cura_SteamEngine %s", VERSION);
			gcode.writeComment("Generated with CubiEngine %s", VERSION);
		}
        return gcode.isOpened();
    }
    
    bool processFile(const char* input_filename)
    {		
#ifdef DEBUG_FORCE_OPTIONS
		//cstyle test
		//config.wipetype = 1;
		//config.wipePosition.X = -70000;//60000;
		//config.wipePosition.Y = 0;//130000;
		//config.wipePosition1.X = 70000;//180000;
		//config.wipePosition1.Y = 0;//130000;
		//config.upSkinCount = 10;
		config.raftType = RAFT_TYPE_CUBICON;// RAFT_TYPE_CUBICON;//RAFT_TYPE_CURA;		
		config.supportInfillType = SUPPORT_TYPE_LINES_END_CONNECT;//SUPPORT_TYPE_GRID;
		config.insetXSpeed = 60;
		config.inset0Speed = 30;
		//config.supportEverywhere = 0;
		config.supportAngle = 60;
		config.skirtLineCount = 0;
		//config.fanSpeedMin = 100;
		//config.fanSpeedMax = 100;
		//config.skirtDistance = 0;
		//config.skirtLineCount = 0;
		//config.skirtMinLength = 0;		
#ifdef DEBUG_CURA_FORCE_OPTIONS
		//config.downSkinCount = 3;
		//config.upSkinCount = 5;
		//config.sparseInfillLineDistance = 3499;
		//config.skirtDistance = 3000;
		//config.skirtLineCount = 3000;
		//config.skirtMinLength = 3000;
		//config.retractionAmountExtruderSwitch = 5000;
		//config.initialLayerSpeed = 25;
		//config.supportLineDistance = 2666;
		//config.minimalLayerTime = 4;
		config.objectPosition.X = 140414;
		config.objectPosition.Y = 113922;
		config.insetEndRetraction = 1;
		//config.retractionZHop = 500;
#endif//DEBUG_CURA_FORCE_OPTIONS
#endif//DEBUG_FORCE_OPTIONS
        if (!gcode.isOpened())
            return false;

        TimeKeeper timeKeeperTotal;
        SliceDataStorage storage;
        preSetup();
        if (!prepareModel(storage, input_filename))
            return false;

		// for FDM skuce
		processSliceData(storage);
		writeGCode(storage);

        logProgress("process", 1, 1);//Report the GUI that a file has been fully processed.
        log("Total time elapsed %5.2fs.\n", timeKeeperTotal.restart());

        return true;
    }
    
    void finalize()
    {
        if (!gcode.isOpened())
            return;
		gcode.finalize(maxObjectHeight, config.moveSpeed, config.endCode, config.finalizeRetraction ? true : false, config.finalizeExtLift);
    }

private:
    void preSetup()
    {        
		skirtConfig.setData(config.skirtSpeed, config.extrusionWidth, "SKIRT");
		homelineConfig.setData(config.skirtSpeed, config.extrusionWidth + config.homelineThickness, "HOMELINE");
        inset0Config.setData(config.inset0Speed, config.extrusionWidth, "WALL-OUTER");
        insetXConfig.setData(config.insetXSpeed, config.extrusionWidth, "WALL-INNER");
        fillConfig.setData(config.infillSpeed, config.extrusionWidth, "FILL");
        supportConfig.setData(config.printSpeed, config.extrusionWidth, "SUPPORT");

        for(unsigned int n=1; n<MAX_EXTRUDERS;n++)
            gcode.setExtruderOffset(n, config.extruderOffset[n].p());
        gcode.setFlavor(config.gcodeFlavor);
        gcode.setRetractionSettings(config.retractionAmount, config.retractionSpeed, config.retractionAmountExtruderSwitch, config.minimalExtrusionBeforeRetraction, config.retractionZHop, config.retractRestoreCorrection);
    }

    bool prepareModel(SliceDataStorage& storage, const char* input_filename)
    {
        timeKeeper.restart();
        log("Loading %s from disk...\n", input_filename);
        SimpleModel* model = NULL;
        if (input_filename[0] == '$')
        {
            //model = new SimpleModel();
            //for(unsigned int n=0; input_filename[n]; n++)
            //{
            //    model->volumes.push_back(SimpleVolume());
            //    SimpleVolume* volume = &model->volumes[model->volumes.size()-1];
            //    //guiSocket.sendNr(GUI_CMD_REQUEST_MESH);
            //    
            //    int32_t vertexCount = guiSocket.recvNr();
            //    int pNr = 0;
            //    log("Reading mesh from socket with %i vertexes\n", vertexCount);
            //    Point3 v[3];
            //    while(vertexCount)
            //    {
            //        float f[3];
            //        guiSocket.recvAll(f, 3 * sizeof(float));
            //        FPoint3 fp(f[0], f[1], f[2]);
            //        v[pNr++] = config.matrix.apply(fp);
            //        if (pNr == 3)
            //        {
            //            volume->addFace(v[0], v[1], v[2]);
            //            pNr = 0;
            //        }
            //        vertexCount--;
            //    }
            //}
        }else{
            model = loadModelFromFile(input_filename, config.matrix);
        }
        if (!model)
        {
            logError("Failed to load model: %s\n", input_filename);
            return false;
        }
        log("Loaded from disk in %5.3fs\n", timeKeeper.restart());
        log("Analyzing and optimizing model...\n");
        OptimizedModel* optimizedModel = new OptimizedModel(model, Point3(config.objectPosition.X, config.objectPosition.Y, -config.objectSink));
        for(unsigned int v = 0; v < model->volumes.size(); v++)
        {
            log("  Face counts: %i -> %i %0.1f%%\n", (int)model->volumes[v].faces.size(), (int)optimizedModel->volumes[v].faces.size(), float(optimizedModel->volumes[v].faces.size()) / float(model->volumes[v].faces.size()) * 100);
            log("  Vertex counts: %i -> %i %0.1f%%\n", (int)model->volumes[v].faces.size() * 3, (int)optimizedModel->volumes[v].points.size(), float(optimizedModel->volumes[v].points.size()) / float(model->volumes[v].faces.size() * 3) * 100);
        }
        delete model;
        log("Optimize model %5.3fs \n", timeKeeper.restart());
        //om->saveDebugSTL("c:\\models\\output.stl");
        
        log("Slicing model...\n");
        vector<Slicer*> slicerList;
        for(unsigned int volumeIdx=0; volumeIdx < optimizedModel->volumes.size(); volumeIdx++)
        {
            //Slicer* slicer = new Slicer(&optimizedModel->volumes[volumeIdx], config.initialLayerThickness - config.layerThickness / 2, config.layerThickness, config.fixHorrible & FIX_HORRIBLE_KEEP_NONE_CLOSED, config.fixHorrible & FIX_HORRIBLE_EXTENSIVE_STITCHING);//org it had a bug related with initail layer thickness
			Slicer* slicer = new Slicer(&optimizedModel->volumes[volumeIdx], config.initialLayerThickness, config.layerThickness, config.fixHorrible & FIX_HORRIBLE_KEEP_NONE_CLOSED, config.fixHorrible & FIX_HORRIBLE_EXTENSIVE_STITCHING);
            slicerList.push_back(slicer);
            for(unsigned int layerNr=0; layerNr<slicer->layers.size(); layerNr++)
            {
                //Reporting the outline here slows down the engine quite a bit, so only do so when debugging.
                //sendPolygonsToGui("outline", layerNr, slicer->layers[layerNr].z, slicer->layers[layerNr].polygonList);
                sendPolygonsToGui("openoutline", layerNr, slicer->layers[layerNr].z, slicer->layers[layerNr].openPolygonList);
            }
        }
        log("Sliced model in %5.3fs\n", timeKeeper.restart());

        log("Generating support map...\n");
        generateSupportGrid(storage.support, optimizedModel, config.supportAngle, config.supportEverywhere > 0, config.supportXYDistance, config.supportZDistance, getRaftThickness());
        
        storage.modelSize = optimizedModel->modelSize;
        storage.modelMin = optimizedModel->vMin;
        storage.modelMax = optimizedModel->vMax;
		Polygon box;
		box.add(Point(storage.modelMin.x, storage.modelMax.y));
		box.add(Point(storage.modelMin.x, storage.modelMin.y));
		box.add(Point(storage.modelMax.x, storage.modelMin.y));
		box.add(Point(storage.modelMax.x, storage.modelMax.y));
		storage.boundaryBox.add(box);		
        delete optimizedModel;

        log("Generating layer parts...\n");
        for(unsigned int volumeIdx=0; volumeIdx < slicerList.size(); volumeIdx++)
        {
            storage.volumes.push_back(SliceVolumeStorage());
            createLayerParts(storage.volumes[volumeIdx], slicerList[volumeIdx], config.fixHorrible & (FIX_HORRIBLE_UNION_ALL_TYPE_A | FIX_HORRIBLE_UNION_ALL_TYPE_B | FIX_HORRIBLE_UNION_ALL_TYPE_C));
            delete slicerList[volumeIdx];
        }
        log("Generated layer parts in %5.3fs\n", timeKeeper.restart());
        return true;
    }
    
    void processSliceData(SliceDataStorage& storage)
    {
        //carveMultipleVolumes(storage.volumes);
        generateMultipleVolumesOverlap(storage.volumes, config.multiVolumeOverlap);
        //dumpLayerparts(storage, "c:/models/output.html");
        
        const unsigned int totalLayers = storage.volumes[0].layers.size();
        for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
        {
            for(unsigned int volumeIdx=0; volumeIdx<storage.volumes.size(); volumeIdx++)
            {
                int insetCount = config.insetCount;
                if (config.spiralizeMode && int(layerNr) < config.downSkinCount && layerNr % 2 == 1)//Add extra insets every 2 layers when spiralizing, this makes bottoms of cups watertight.
                    insetCount += 5;
                SliceLayer* layer = &storage.volumes[volumeIdx].layers[layerNr];
                generateInsets(layer, config.extrusionWidth, insetCount);

                for(unsigned int partNr=0; partNr<layer->parts.size(); partNr++)
                {
                    if (layer->parts[partNr].insets.size() > 0)
                    {
                        sendPolygonsToGui("inset0", layerNr, layer->z, layer->parts[partNr].insets[0]);
                        for(unsigned int inset=1; inset<layer->parts[partNr].insets.size(); inset++)
                            sendPolygonsToGui("insetx", layerNr, layer->z, layer->parts[partNr].insets[inset]);
                    }
                }
            }
            logProgress("inset",layerNr+1,totalLayers);
        }
        if (config.enableOozeShield)
        {
            for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
            {
                Polygons oozeShield;
                for(unsigned int volumeIdx=0; volumeIdx<storage.volumes.size(); volumeIdx++)
                {
                    for(unsigned int partNr=0; partNr<storage.volumes[volumeIdx].layers[layerNr].parts.size(); partNr++)
                    {
                        oozeShield = oozeShield.unionPolygons(storage.volumes[volumeIdx].layers[layerNr].parts[partNr].outline.offset(2000));
                    }
                }
                storage.oozeShield.push_back(oozeShield);
            }
            
            for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
                storage.oozeShield[layerNr] = storage.oozeShield[layerNr].offset(-1000).offset(1000);
            int offsetAngle = tan(60.0*M_PI/180) * config.layerThickness;//Allow for a 60deg angle in the oozeShield.
            for(unsigned int layerNr=1; layerNr<totalLayers; layerNr++)
                storage.oozeShield[layerNr] = storage.oozeShield[layerNr].unionPolygons(storage.oozeShield[layerNr-1].offset(-offsetAngle));
            for(unsigned int layerNr=totalLayers-1; layerNr>0; layerNr--)
                storage.oozeShield[layerNr-1] = storage.oozeShield[layerNr-1].unionPolygons(storage.oozeShield[layerNr].offset(-offsetAngle));
        }
        log("Generated inset in %5.3fs\n", timeKeeper.restart());
		
        for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
        {
            if (!config.spiralizeMode || int(layerNr) < config.downSkinCount)    //Only generate up/downskin and infill for the first X layers when spiralize is choosen.
            {
                for(unsigned int volumeIdx=0; volumeIdx<storage.volumes.size(); volumeIdx++)
                {					
					int topMarginCount = config.upSkinCount>0 ? config.layerSupportTopMarginCount:0;
					if (config.layerSupport || layerNr<=config.downSkinCount || layerNr>=(totalLayers-(config.upSkinCount+topMarginCount)))//cstyle bridge? test. generate only up/down skins not middle skins.
					generateSkins(layerNr, storage.volumes[volumeIdx], config.extrusionWidth, config.downSkinCount, config.upSkinCount, config.infillOverlap);
                    generateSparse(layerNr, storage.volumes[volumeIdx], config.extrusionWidth, config.downSkinCount, config.upSkinCount);
                    
                    SliceLayer* layer = &storage.volumes[volumeIdx].layers[layerNr];
                    for(unsigned int partNr=0; partNr<layer->parts.size(); partNr++)
                        sendPolygonsToGui("skin", layerNr, layer->z, layer->parts[partNr].skinOutline);
                }
            }
            logProgress("skin",layerNr+1,totalLayers);
        }
        log("Generated up/down skin in %5.3fs\n", timeKeeper.restart());

		//cstyle modified. add cleaning head code.
		//if (config.wipeTowerSize > 0)
        if (config.wipetype == 0 && config.wipeTowerSize > 0)
        {
            PolygonRef p = storage.wipeTower.newPoly();
            p.add(Point(storage.modelMin.x - 3000, storage.modelMax.y + 3000));
            p.add(Point(storage.modelMin.x - 3000, storage.modelMax.y + 3000 + config.wipeTowerSize));
            p.add(Point(storage.modelMin.x - 3000 - config.wipeTowerSize, storage.modelMax.y + 3000 + config.wipeTowerSize));
            p.add(Point(storage.modelMin.x - 3000 - config.wipeTowerSize, storage.modelMax.y + 3000));
            
            storage.wipePoint = Point(storage.modelMin.x - 3000 - config.wipeTowerSize / 2, storage.modelMax.y + 3000 + config.wipeTowerSize / 2);
        } else
		if (config.wipetype==1)
		{
			PolygonRef p = storage.wipeTower.newPoly();//not used here
			storage.wipePoint = Point(config.wipePosition.X, config.wipePosition.Y);
			storage.wipePoint1 = Point(config.wipePosition1.X, config.wipePosition1.Y);
			//cstyle test
			/*Point center = Point((storage.modelMin.x +storage.modelMax.x) / 2, (storage.modelMax.y + storage.modelMin.y) / 2);
			Point offset(80000, 0);			
			storage.wipePoint = center-offset;			
			storage.wipePoint1 = center+offset;*/
		}
				
        generateSkirt(storage, config.skirtDistance, config.extrusionWidth, config.skirtLineCount, config.skirtMinLength, config.initialLayerThickness);
		generateRaft(storage, config.raftMargin, config.raftDashStride, config.raftDashSize);
		sendPolygonsToGui("skirt", 0, config.initialLayerThickness, storage.skirt);
        
        for(unsigned int volumeIdx=0; volumeIdx<storage.volumes.size(); volumeIdx++)
        {
            for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
            {
                for(unsigned int partNr=0; partNr<storage.volumes[volumeIdx].layers[layerNr].parts.size(); partNr++)
                {
                    if (layerNr > 0)
                        storage.volumes[volumeIdx].layers[layerNr].parts[partNr].bridgeAngle = bridgeAngle(&storage.volumes[volumeIdx].layers[layerNr].parts[partNr], &storage.volumes[volumeIdx].layers[layerNr-1]);
                    else
                        storage.volumes[volumeIdx].layers[layerNr].parts[partNr].bridgeAngle = -1;
                }
            }
        }
    }
	
    void writeGCode(SliceDataStorage& storage)
    {
        if (fileNr == 1)
        {
            if (gcode.getFlavor() == GCODE_FLAVOR_ULTIGCODE)
            {
                gcode.writeCode(";FLAVOR:UltiGCode");
                gcode.writeCode(";TIME:<__TIME__>");
                gcode.writeCode(";MATERIAL:<FILAMENT>");
                gcode.writeCode(";MATERIAL2:<FILAMEN2>");
            }
            gcode.writeCode(config.startCode);
        } else
		{
            gcode.writeFanCommand(0);
            gcode.resetExtrusionValue();
            gcode.writeRetraction();
            gcode.setZ(maxObjectHeight + 5000);
            gcode.writeMove(Point(storage.modelMin.x, storage.modelMin.y), config.moveSpeed, 0);
        }	

		//cstyle
		//writeOptions(gcode, config);		
		
        fileNr++;
        
		//cstyle
		int raftThickness = 0;
		int nLayerMiddle = 2;
		int airGap = 0;//it is defined as initialLayerThickness
		bool useAirGap = false;
		if (storage.volumes[0].layers[0].parts.size() == 0)
			log("error\fzero_layer\n");
		storage.volumes[0].RemoveEmptyLayer();
        unsigned int totalLayers = storage.volumes[0].layers.size();
		

        gcode.writeComment("Layer count: %d", totalLayers);
        if (config.raftBaseThickness > 0 && config.raftInterfaceThickness > 0)
        {
			GCodePathConfig raftBaseConfig(config.initialLayerSpeed, config.raftBaseLinewidth, "RAFT");
            GCodePathConfig raftInterfaceConfig(config.initialLayerSpeed, config.raftInterfaceLinewidth, "RAFT");

			if (config.raftType==RAFT_TYPE_CURA)
			{
				{
					gcode.writeComment("LAYER:-2");
					gcode.writeComment("TYPE:RAFT");
					GCodePlanner gcodeLayer(gcode, config.moveSpeed, config.retractionMinimalDistance);
					gcode.setZ(config.raftBaseThickness);
					gcode.setExtrusion(config.raftBaseThickness, config.filamentDiameter, config.filamentFlow);
					gcodeLayer.addPolygonsByOptimizer(storage.raftOutline, &raftBaseConfig);
                
					Polygons raftLines;
					generateLineInfill(storage.raftOutline, raftLines, config.raftBaseLinewidth, config.raftLineSpacing, config.infillOverlap, 0);
					gcodeLayer.addPolygonsByOptimizer(raftLines, &raftBaseConfig);
                
					gcodeLayer.writeGCode(false, config.raftBaseThickness);
				}
				{
					gcode.writeComment("LAYER:-1");
					gcode.writeComment("TYPE:RAFT");
					GCodePlanner gcodeLayer(gcode, config.moveSpeed, config.retractionMinimalDistance);
					gcode.setZ(config.raftBaseThickness + config.raftInterfaceThickness);
					gcode.setExtrusion(config.raftInterfaceThickness, config.filamentDiameter, config.filamentFlow);
                
					Polygons raftLines;
					generateLineInfill(storage.raftOutline, raftLines, config.raftInterfaceLinewidth, config.raftLineSpacing, config.infillOverlap, 90);
					gcodeLayer.addPolygonsByOptimizer(raftLines, &raftInterfaceConfig);
                
					gcodeLayer.writeGCode(false, config.raftInterfaceThickness);
				}
			} else
			if (config.raftType==RAFT_TYPE_CUBICON)
			{
#ifdef DEBUG_FORCE_OPTIONS
				//test force value			
				//config.raftBaseThickness = 500;
				//config.raftBaseLinewidth = 700;
				//config.raftInterfaceThickness = 300;
				//config.raftInterfaceLinewidth = 300;
#endif	
				//speed
				raftBaseConfig.speed = config.raftBaseSpeed;
				//config.raftInterfaceSpeed = config.raftBaseSpeed;
				raftInterfaceConfig.speed = config.raftInterfaceSpeed;

				char str[128];
				int layerNr=0;
				bool useBaseRaft = true;
				bool useDiagonalRaft = true;
				bool useMiddleRaft = true;
				useAirGap = true;

				//adjust temperature of raft for TPU filament.
				LayerControl layerCtrl;
				if (GetLayerControl(-5, &layerCtrl))
					gcode.writeLine(layerCtrl.gcode.c_str());				

				//base raft
				if (useBaseRaft)
				{
					layerNr = 1 + (useDiagonalRaft?1:0) + (useMiddleRaft?nLayerMiddle:0) +  (useAirGap?1:0);
					raftThickness = writeRaft(gcode, config, storage, raftBaseConfig, layerNr, 
						PATTERN_ZIGZAG, 0, raftThickness, config.raftBaseThickness,						
						config.raftBaseLinewidth, config.raftLineSpacing, true, 1500, true, true);
				}
				//base diagonal raft
				if (useDiagonalRaft)
				{												
					config.raftLineSpacing = config.raftInterfaceLinewidth*4;
					layerNr = 1 + (useMiddleRaft?nLayerMiddle:0) +  (useAirGap?1:0);
					raftThickness = writeRaft(gcode, config, storage, raftInterfaceConfig, layerNr, 
						PATTERN_LINE, 45, raftThickness, config.raftInterfaceThickness, 
						config.raftInterfaceLinewidth, config.raftLineSpacing, false, 0);
				}
				//middle raft
				if (useMiddleRaft)
				{
					int lineThickness = 0;
					int lineWidth = 0;
					int lineSpacing = 0;
					for (int i=0;i<nLayerMiddle;i++)
					{
						if (i==0)
						{						
							lineThickness = config.raftPlaneThickness1;
							lineWidth = config.raftPlaneWidth1;
							lineSpacing = config.raftInterfaceLinewidth;
						} if (i==1)
						{	
							lineThickness = config.raftPlaneThickness2;
							lineWidth = config.raftPlaneWidth2;
							lineSpacing = config.raftInterfaceLinewidth;
						}
						int layerNr = nLayerMiddle-i + (useAirGap?1:0);
						raftThickness = writeRaft(gcode, config, storage, raftInterfaceConfig, layerNr, 
							PATTERN_LINE, (i % 2 == 0) ? 90 : 0, raftThickness, lineThickness,
						lineWidth, lineSpacing, false, 0);
					}
				}				
			}
        }

        int volumeIdx = 0;
		Point startP = gcode.getPositionXY();
        for(unsigned int layerNr=0; layerNr<totalLayers; layerNr++)
        {
            logProgress("export", layerNr+1, totalLayers);

            if (int(layerNr) < config.initialSpeedupLayers)
            {
                int n = config.initialSpeedupLayers;
                //skirtConfig.setData(config.printSpeed * layerNr / n + config.initialLayerSpeed * (n - layerNr) / n, config.extrusionWidth, "SKIRT");
				skirtConfig.setData(config.skirtSpeed, config.extrusionWidth, "SKIRT");
                inset0Config.setData(config.inset0Speed * layerNr / n + config.initialLayerSpeed * (n - layerNr) / n, config.extrusionWidth, "WALL-OUTER");
                insetXConfig.setData(config.insetXSpeed * layerNr / n + config.initialLayerSpeed * (n - layerNr) / n, config.extrusionWidth, "WALL-INNER");
                fillConfig.setData(config.infillSpeed * layerNr / n + config.initialLayerSpeed * (n - layerNr) / n, config.extrusionWidth, "FILL");
                supportConfig.setData(config.printSpeed * layerNr / n + config.initialLayerSpeed * (n - layerNr) / n, config.extrusionWidth, "SUPPORT");
            }else{
                //skirtConfig.setData(config.printSpeed, config.extrusionWidth, "SKIRT");
				skirtConfig.setData(config.skirtSpeed, config.extrusionWidth, "SKIRT");
                inset0Config.setData(config.inset0Speed, config.extrusionWidth, "WALL-OUTER");
                insetXConfig.setData(config.insetXSpeed, config.extrusionWidth, "WALL-INNER");
                fillConfig.setData(config.infillSpeed, config.extrusionWidth, "FILL");
                supportConfig.setData(config.printSpeed, config.extrusionWidth, "SUPPORT");
            }
			
            gcode.writeComment("LAYER:%d", layerNr);
			//cstyle layer control
			LayerControl layerCtrl;
			if (GetLayerControl(layerNr, &layerCtrl))
			{
				gcode.writeLine(layerCtrl.gcode.c_str());
			}

            if (layerNr == 0)
                gcode.setExtrusion(config.initialLayerThickness, config.filamentDiameter, config.filamentFlow);
            else
                gcode.setExtrusion(config.layerThickness, config.filamentDiameter, config.filamentFlow);

            //GCodePlanner gcodeLayer(gcode, config.moveSpeed, config.retractionMinimalDistance);//org
			//set force the first travel speed.		
			GCodePlanner gcodeLayer(gcode, config.moveSpeed, config.retractionMinimalDistance, (layerNr==0  && !isUsedRaft())?config.startSpeed:0);//cstyle

			//cstyle airgap
			int airGap = 0;
			if (layerNr == 0 && useAirGap && config.airGap > 0)
			{
				airGap = config.airGap;
			}

            int32_t z = config.initialLayerThickness + layerNr * config.layerThickness + raftThickness + airGap;            
            gcode.setZ(z);			
            
            bool printSupportFirst = (storage.support.generated && config.supportExtruder > 0 && config.supportExtruder == gcodeLayer.getExtruder());
            if (printSupportFirst)
                addSupportToGCode(storage, gcodeLayer, layerNr, raftThickness);
            

			for(unsigned int volumeCnt = 0; volumeCnt < storage.volumes.size(); volumeCnt++)
            {				
				if (volumeCnt > 0)
					volumeIdx = (volumeIdx + 1) % storage.volumes.size();
				addVolumeLayerToGCode(storage, gcodeLayer, volumeIdx, layerNr, NULL);//, &startP);				
            }
            if (!printSupportFirst)
                addSupportToGCode(storage, gcodeLayer, layerNr, raftThickness);
            
            //Finish the layer by applying speed corrections for minimal layer times
            gcodeLayer.forceMinimalLayerTime(config.minimalLayerTime, config.minimalFeedrate);

            int fanSpeed = config.fanSpeedMin;
            if (gcodeLayer.getExtrudeSpeedFactor() <= 50)
            {
                fanSpeed = config.fanSpeedMax;
            }else{
                int n = gcodeLayer.getExtrudeSpeedFactor() - 50;
                fanSpeed = config.fanSpeedMin * n / 50 + config.fanSpeedMax * (50 - n) / 50;
            }
            if (int(layerNr) < config.fanFullOnLayerNr)
            {
                //Slow down the fan on the layers below the [fanFullOnLayerNr], where layer 0 is speed 0.
                fanSpeed = fanSpeed * layerNr / config.fanFullOnLayerNr;
            }
			//cstyle. there are many kind of fan motor.  generally it is controlled by pwm motor controller. 
			//but a some motor would not work at speed 1 or low speed, because of it's own threshold voltage which enables the motor to rotate.
			//To solve this problem, I applied the fanspeedThreshold in the settingConfig.
			if (config.fanspeedThreshold>0 && fanSpeed>0)			
				fanSpeed = readjustFanspeed(fanSpeed, config.fanspeedThreshold);
			
            gcode.writeFanCommand(fanSpeed);
			gcodeLayer.writeGCode(config.coolHeadLift > 0, int(layerNr) > 0 ? config.layerThickness : config.initialLayerThickness, (layerNr == 0) ? true : false);
        }
        
        log("Wrote layers in %5.2fs.\n", timeKeeper.restart());
        gcode.tellFileSize();
        gcode.writeFanCommand(0);

        //Store the object height for when we are printing multiple objects, as we need to clear every one of them when moving to the next position.
		//cstyle
        maxObjectHeight = std::max(maxObjectHeight, storage.modelSize.z - config.objectSink);
    }
		   
    //Add a single layer from a single mesh-volume to the GCode
    void addVolumeLayerToGCode(SliceDataStorage& storage, GCodePlanner& gcodeLayer, int volumeIdx, int layerNr, Point *pStartPoint)
    {
        int prevExtruder = gcodeLayer.getExtruder();
        bool extruderChanged = gcodeLayer.setExtruder(volumeIdx);		
		if (!isUsedRaft() && layerNr == 0 && volumeIdx == 0 && storage.skirt.size())
		{			
			addHomeLine(gcodeLayer, storage.skirt);//cstyle			
			gcodeLayer.addPolygonsByOptimizer(storage.skirt, &skirtConfig);
		}

        SliceLayer* layer = &storage.volumes[volumeIdx].layers[layerNr];
		//org
		//if (extruderChanged)
		//	addWipeTower(storage, gcodeLayer, layerNr, prevExtruder);
		if (extruderChanged && config.wipetype==0) 				
			addWipeTower(storage, gcodeLayer, layerNr, prevExtruder);
		else
		if (extruderChanged && config.wipetype==1 && layer && layer->parts.size())		
			addWipeExtruder(storage, gcodeLayer, layerNr, prevExtruder, (double)config.wipeExtrusionAmount);		
        
        if (storage.oozeShield.size() > 0 && storage.volumes.size() > 1)
        {
			if (layerNr == 0)			
				addHomeLine(gcodeLayer, storage.oozeShield[layerNr]);//cstyle	
			
            gcodeLayer.setAlwaysRetract(true);
            gcodeLayer.addPolygonsByOptimizer(storage.oozeShield[layerNr], &skirtConfig);
            sendPolygonsToGui("oozeshield", layerNr, layer->z, storage.oozeShield[layerNr]);
            gcodeLayer.setAlwaysRetract(!config.enableCombing);
        }
        
        PathOrderOptimizer partOrderOptimizer(gcode.getPositionXY());
        for(unsigned int partNr=0; partNr<layer->parts.size(); partNr++)
        {
            partOrderOptimizer.addPolygon(layer->parts[partNr].insets[0][0]);
        }
        partOrderOptimizer.optimize();
        
        for(unsigned int partCounter=0; partCounter<partOrderOptimizer.polyOrder.size(); partCounter++)
        {
            SliceLayerPart* part = &layer->parts[partOrderOptimizer.polyOrder[partCounter]];
            
            if (config.enableCombing)
                gcodeLayer.setCombBoundary(&part->combBoundery);
            else
                gcodeLayer.setAlwaysRetract(true);
            
            if (config.insetCount > 0)
            {
                if (config.spiralizeMode)
                {
                    if (int(layerNr) >= config.downSkinCount)
                        inset0Config.spiralize = true;
					if (int(layerNr) == config.downSkinCount && part->insets.size() > 0)
					{
						if (layerNr == 0)
							addHomeLine(gcodeLayer, part->insets[0]);//cstyle	
						gcodeLayer.addPolygonsByOptimizer(part->insets[0], &insetXConfig);
					}
                }

				//cstyle it decides the order that drawing a wall.
				//2016.05.03 changed that retraction is occurring only outline wall.
				//2016.06.03 applied gcode path control(endReduce, limit endRetraction, endGap, speedMatching, forceEndPath) and merged in->out and out->in code to one for loop.
				bool retract = config.insetEndRetraction?true:false;
				Point lastPoint;
				
				for (int i = 0; i < part->insets.size();i++)
				{
					int insetNr = i;//outline -> inline
					if (config.wallInOutOrder == 0)//incompositeline -> outline	 -- default value
						insetNr = part->insets.size() - 1 - i;
					//cstyle modified to start almost the same position . 
					//but actually each polygon's start point are different, so it's impossible to have the same start point for all layers.
					if (i == 0)
						addHomeLine(gcodeLayer, part->insets[insetNr]);//cstyle
					
					if (config.fixedWallBeginAngle!=-999)
						RearrangeStartPosition(part->insets[insetNr], config.fixedWallBeginAngle);

					if (insetNr == 0)
					{							
						GCodeControl ctrl;
						ctrl.endRetraction = false;
						ctrl.travelSpeedMatching = config.travelSpeedMatching==1?true:false;	 
						ctrl.endGap = config.outerWallEndGap;
						ctrl.outwallStartOffset = config.outerWallStartOffset;
						if (config.fixedWallBeginAngle != -999)
							ctrl.useSameStartPosition = true;
																		
						gcodeLayer.addPolygonsByOptimizer(part->insets[insetNr], &inset0Config, ctrl);//outer wall

						//cstyle pull end of wall by travel
						if (config.pullEndOfWallTravelDistance > 0)
						{
							GCodePath *pPath = gcodeLayer.getLatestPathWithConfig(&inset0Config);
							int size = pPath->points.size();
							if (size > 1)
							{
								Vector2 v2 = FindPerpendicularVector(pPath->points[size - 1], pPath->points[size - 2], DEGTORAD((float)config.pullEndOfWallAngle), config.pullEndOfWallTravelDistance);
								Point perpendicularP = Point((int)v2.X, (int)v2.Y);
								gcodeLayer.forceAddTravel(perpendicularP, ctrl.travelSpeedMatching ? inset0Config.speed : 0);
							}
						}
					}
					else
					{
						GCodeControl ctrl;
						if (config.fixedWallBeginAngle != -999)
							ctrl.useSameStartPosition = true;
						ctrl.endRetraction = retract;
						if (config.innerWallEndControlType==0)
							ctrl.endGap = (insetNr == 1)?config.innerWallEndGap:0;
						else
						{
							if (insetNr == 1)
							{
								ctrl.endReduce = true;
								ctrl.reduceDist = config.innerWallEndGap;
								ctrl.reduceRate = (float)config.innerwallEndReduceRate/100.0f;
							}
							else
								ctrl.endReduce = false;
						}
						gcodeLayer.addPolygonsByOptimizer(part->insets[insetNr], &insetXConfig, ctrl);//inner_wall
					}
				}
            }
            
            Polygons fillPolygons;
            int fillAngle = 45;
            if (layerNr & 1) 
                fillAngle += 90;

			if (config.initialLayerFillType == 1 && layerNr == 0)
			{
				//fillAngle = 90;
				int sparseSteps[1] = { config.extrusionWidth };
				generateConcentricInfill(part->skinOutline, fillPolygons, sparseSteps, 1);
			}
			else
            //cstyle - bridge layerFill?	
				generateLineInfill(part->skinOutline, fillPolygons, config.extrusionWidth, config.extrusionWidth, config.infillOverlap, (part->bridgeAngle > -1) ? part->bridgeAngle : fillAngle);

            //int sparseSteps[2] = {config.extrusionWidth*5, config.extrusionWidth * 0.8};
            //generateConcentricInfill(part->sparseOutline, fillPolygons, sparseSteps, 2);
						
            if (config.sparseInfillLineDistance > 0)
            {
				if (config.sparseInfillLineDistance > config.extrusionWidth * 4)
                {
                    generateLineInfill(part->sparseOutline, fillPolygons, config.extrusionWidth, config.sparseInfillLineDistance * 2, config.infillOverlap, 45);
                    generateLineInfill(part->sparseOutline, fillPolygons, config.extrusionWidth, config.sparseInfillLineDistance * 2, config.infillOverlap, 45 + 90);
                } else
                {
                    generateLineInfill(part->sparseOutline, fillPolygons, config.extrusionWidth, config.sparseInfillLineDistance, config.infillOverlap, fillAngle);
                }
            }

            gcodeLayer.addPolygonsByOptimizer(fillPolygons, &fillConfig);
            sendPolygonsToGui("infill", layerNr, layer->z, fillPolygons);
            
            //After a layer part, make sure the nozzle is inside the comb boundary, so we do not retract on the perimeter.
            if (!config.spiralizeMode || int(layerNr) < config.downSkinCount)
                gcodeLayer.moveInsideCombBoundary(config.extrusionWidth * 2);
        }
        gcodeLayer.setCombBoundary(NULL);
    }
    
    void addSupportToGCode(SliceDataStorage& storage, GCodePlanner& gcodeLayer, int layerNr, int raftThickness)
    {
        if (!storage.support.generated)
            return;
        
        if (config.supportExtruder > -1)
        {
            int prevExtruder = gcodeLayer.getExtruder();
            if (gcodeLayer.setExtruder(config.supportExtruder))
			{                
				if (config.wipetype==0)
					addWipeTower(storage, gcodeLayer, layerNr, prevExtruder);
				else
					addWipeExtruder(storage, gcodeLayer, layerNr, prevExtruder, (double)config.wipeExtrusionAmount);
			}
            
            if (storage.oozeShield.size() > 0 && storage.volumes.size() == 1)
            {
                gcodeLayer.setAlwaysRetract(true);
                gcodeLayer.addPolygonsByOptimizer(storage.oozeShield[layerNr], &skirtConfig);
                gcodeLayer.setAlwaysRetract(!config.enableCombing);
            }
        }
        int32_t z = config.initialLayerThickness + layerNr * config.layerThickness;
		z += raftThickness;
        SupportPolyGenerator supportGenerator(storage.support, z);
        for(unsigned int volumeCnt = 0; volumeCnt < storage.volumes.size(); volumeCnt++)
        {
            SliceLayer* layer = &storage.volumes[volumeCnt].layers[layerNr];
            for(unsigned int n=0; n<layer->parts.size(); n++)
                supportGenerator.polygons = supportGenerator.polygons.difference(layer->parts[n].outline.offset(config.supportXYDistance));
        }
        //Contract and expand the suppory polygons so small sections are removed and the final polygon is smoothed a bit.
        supportGenerator.polygons = supportGenerator.polygons.offset(-config.extrusionWidth * 3);
        supportGenerator.polygons = supportGenerator.polygons.offset(config.extrusionWidth * 3);
        sendPolygonsToGui("support", layerNr, z, supportGenerator.polygons);
        
        vector<Polygons> supportIslands = supportGenerator.polygons.splitIntoParts();
        
        PathOrderOptimizer islandOrderOptimizer(gcode.getPositionXY());
        for(unsigned int n=0; n<supportIslands.size(); n++)
        {
            islandOrderOptimizer.addPolygon(supportIslands[n][0]);
        }
        islandOrderOptimizer.optimize();
        
        for(unsigned int n=0; n<supportIslands.size(); n++)
        {
            Polygons& island = supportIslands[islandOrderOptimizer.polyOrder[n]];
            
            Polygons supportLines;
            if (config.supportLineDistance > 0)
            {
				if (config.supportInfillType==SUPPORT_TYPE_GRID)
				{
					if (config.supportLineDistance > config.extrusionWidth * 4)
					{
						generateLineInfill(island, supportLines, config.extrusionWidth, config.supportLineDistance*2, config.infillOverlap, 0);
						generateLineInfill(island, supportLines, config.extrusionWidth, config.supportLineDistance*2, config.infillOverlap, 90);
					}else
					{
						generateLineInfill(island, supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap, (layerNr & 1) ? 0 : 90);
					}
				} else
				if (config.supportInfillType==SUPPORT_TYPE_LINES)
				{
					generateLineInfillEx(supportIslands[n], supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap, 
						0, PATTERN_LINE, 0, 0);
				} else
				if (config.supportInfillType==SUPPORT_TYPE_LINES_END_CONNECT)
				{
					Polygons refInsets;
					for (int i=0;i<storage.volumes.size();i++)
					{
						SliceLayer* layer = &storage.volumes[i].layers[layerNr];
						for(unsigned int partNr=0; partNr<layer->parts.size(); partNr++)
						{
							if (layer->parts[partNr].insets.size() > 0)
								refInsets.add(layer->parts[partNr].insets[0].offset(config.extrusionWidth));
						}
					}

					Polygons *pBoundary = &storage.boundaryBox;	
					if (layerNr < config.strongSupportLayerNr)//grid type. this makes strong base
					{
						generateLineInfillEx(island, supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap,
							0, PATTERN_LINE, 0, 0, pBoundary, (refInsets.size() > 0) ? &refInsets : NULL);
						generateLineInfillEx(island, supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap,
							90, PATTERN_LINE, 0, 0, pBoundary, (refInsets.size() > 0) ? &refInsets : NULL);						
					}
					else
					{						
						//generateLineInfillEx(supportIslands[n], supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap,
						if (layerNr%30<2)
							generateLineInfillEx(island, supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap,
							0, PATTERN_LINE, 1, 1, pBoundary, (refInsets.size() > 0) ? &refInsets : NULL);
						else
							generateLineInfillEx(island, supportLines, config.extrusionWidth, config.supportLineDistance, config.infillOverlap,
							0, PATTERN_LINE, 1, 0, pBoundary, (refInsets.size() > 0) ? &refInsets : NULL);
					}
				}
			}

			//cstyle test
			//if (storage.volumes.size() > 0 && layerNr<storage.volumes[0].layers.size()-3)
			//{
			//	Polygons *pBoundary = &storage.boundaryBox;
			//	SliceLayer* nextLayer = &storage.volumes[0].layers[layerNr+3];
			//	SliceLayer* layer = &storage.volumes[0].layers[layerNr];
			//	if (nextLayer != NULL && layer != NULL && nextLayer->parts.size()>0 && layer->parts.size()>0)
			//	{
			//		Polygons diff = nextLayer->parts[0].outline.difference(layer->parts[0].outline);
			//		//Polygons inter = island.intersection(diff);
			//		Polygons base = diff.offset(-600);
			//		Polygons diffInfill;
			//		for (int i = 0; i < base.size(); )
			//		{
			//			if (fabs(base[i].area()) < 5000000.0f)
			//			{
			//				base.remove(i);
			//				i = 0;
			//			} else
			//			  i++;
			//		}					
			//		generateLineInfillEx(base, diffInfill, config.extrusionWidth, config.supportLineDistance / 5, config.infillOverlap,
			//			90, PATTERN_LINE, 0, 0, pBoundary, NULL);// (refInsets.size() > 0) ? &refInsets : NULL);
			//		gcodeLayer.addPolygonsByOptimizer(diffInfill, &supportConfig);
			//		gcodeLayer.addPolygonsByOptimizer(base, &supportConfig);
			//	}
			//}
        
            gcodeLayer.forceRetract();
            if (config.enableCombing)
                gcodeLayer.setCombBoundary(&island);
			
			//int strongSupportLayerNr = 2;
			//int strongLayerCnt = layerNr % 30;
			//if (strongLayerCnt<strongSupportLayerNr || config.supportInfillType == SUPPORT_TYPE_GRID)
			if (layerNr<2 || config.supportInfillType == SUPPORT_TYPE_GRID)
				gcodeLayer.addPolygonsByOptimizer(island, &supportConfig);//org

            gcodeLayer.addPolygonsByOptimizer(supportLines, &supportConfig);
            gcodeLayer.setCombBoundary(NULL);
        }
    }
    
    void addWipeTower(SliceDataStorage& storage, GCodePlanner& gcodeLayer, int layerNr, int prevExtruder)
    {
        if (config.wipeTowerSize < 1)
            return;
        //If we changed extruder, print the wipe/prime tower for this nozzle;
        gcodeLayer.addPolygonsByOptimizer(storage.wipeTower, &supportConfig);
        Polygons fillPolygons;
        generateLineInfill(storage.wipeTower, fillPolygons, config.extrusionWidth, config.extrusionWidth, config.infillOverlap, 45 + 90 * (layerNr % 2));
        gcodeLayer.addPolygonsByOptimizer(fillPolygons, &supportConfig);
        
        //Make sure we wipe the old extruder on the wipe tower.
        gcodeLayer.addTravel(storage.wipePoint - config.extruderOffset[prevExtruder].p() + config.extruderOffset[gcodeLayer.getExtruder()].p());
    }

	//cstyle additional functions
	void addWipeExtruder(SliceDataStorage& storage, GCodePlanner& gcodeLayer, int layerNr, int prevExtruder, double extrusionAmount)
    {
		Point wipePoint = storage.wipePoint;
		if (gcodeLayer.getExtruder()==1)
			wipePoint = storage.wipePoint1;
		gcodeLayer.addWipeExtrusion(wipePoint, extrusionAmount);
        //Make sure we wipe the old extruder on the wipe tower.
        //gcodeLayer.addTravel(storage.wipePoint - config.extruderOffset[prevExtruder].p() + config.extruderOffset[gcodeLayer.getExtruder()].p());
		//gcodeLayer.addTravel(wipePoint - config.extruderOffset[prevExtruder].p() + config.extruderOffset[gcodeLayer.getExtruder()].p());
    }

	void changeExtruderWithWipeTower(GCodePlanner &gcodeLayer, int extruderNr, SliceDataStorage &storage, int layerNr)
	{
		int prevExtruder = gcodeLayer.getExtruder();
		if (gcodeLayer.setExtruder(extruderNr))
		{
			if (config.wipetype==0)
				addWipeTower(storage, gcodeLayer, layerNr, prevExtruder);
			else
			if (config.wipetype==1)
				addWipeExtruder(storage, gcodeLayer, layerNr, prevExtruder, config.wipeExtrusionAmount);
		}
	}

	bool isUsedRaft()
	{
		if (config.raftBaseThickness > 0 && config.raftInterfaceThickness > 0)
			return true;
		return false;
	}

	int getRaftThickness()
	{
		if (isUsedRaft())
			return config.raftBaseThickness + config.raftInterfaceThickness + config.raftPlaneThickness1 + config.raftPlaneThickness2;
		return 0;
	}

	int writeRaft(GCodeExport &gcode, 
		ConfigSettings &config, 
		SliceDataStorage& storage, 
		GCodePathConfig &pathConfig, 
		int layerNr, 
		PATTERN_TYPE patternType, 
		double rotation, 
		int raftThickness, 
		int thickness, 
		int lineWidth, 
		int lineSpacing, 
		bool outline, 
		int arg, 
		bool optimize=false,
		bool useHomeLine=false)
	{
		char str[1024];
		sprintf(str, "LAYER:-%d\0", layerNr);
		gcode.writeComment(str);
		//gcode.writeComment("TYPE:RAFT");
				
		GCodePlanner gcodeLayer(gcode, config.moveSpeed, config.retractionMinimalDistance, (layerNr==5)?config.startSpeed:0);
		//cstyle changes the support extruder if it supporter extruder is assigned.
		if (config.supportExtruder > -1)
			changeExtruderWithWipeTower(gcodeLayer, config.supportExtruder, storage, layerNr);
		
		raftThickness += thickness;
		gcode.setZ(raftThickness);	
		gcode.setExtrusion(thickness, config.filamentDiameter, config.filamentFlow);
		
		if (outline)
		{
			if (useHomeLine)
				addHomeLine(gcodeLayer, storage.raftOutline);
			gcodeLayer.addPolygonsByOptimizer(storage.raftOutline, &pathConfig);
		}

		//add raft infill
		Polygons raftInfill;
		generateLineInfillEx(storage.raftOutline, raftInfill, lineWidth, lineSpacing, config.infillOverlap, rotation, patternType, arg);
		if (optimize)
		{
			Polygons optimizedInfill;
			optimizePolylinePath(optimizedInfill, raftInfill, Point(0,0));
			GCodeControl ctrl;
			ctrl.endRetraction = true;
			ctrl.keepPointOrder = true;
			gcodeLayer.addPolygonsByOptimizer(optimizedInfill, &pathConfig, ctrl);// true, true);
		} else
			gcodeLayer.addPolygonsByOptimizer(raftInfill, &pathConfig);
		
		gcodeLayer.writeGCode(false, config.raftBaseThickness);
		return raftThickness;
	}

	//cstyle
	void addHomeLine(GCodePlanner &gcodeLayer, Polygons &polygons)
	{
		if (drawHomeLine)
		{
			Polygons firstLine = genInitPath(config, 50000, Point(0, 0), &polygons);
			GCodeControl ctrl;
			ctrl.endRetraction = false;
			ctrl.keepPointOrder = true;
			gcodeLayer.addPolygonsByOptimizer(firstLine, &homelineConfig, ctrl);// false, true);
			drawHomeLine = false;
		}
	}

	bool GetLayerControl(int layerIndex, LayerControl *outData)
	{
		std::vector<LayerControl> &layerControls = config.layerControls;
		for (int i = 0; i < layerControls.size(); i++)
		{			
			LayerControl data = layerControls[i];
			if (data.layerIndex != -1 && data.layerIndex == layerIndex)
			{
				*outData = data;
				return true;
			}
		}
		return false;
	}

};

#endif//FFF_PROCESSOR_H
