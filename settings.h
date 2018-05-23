#ifndef SETTINGS_H
#define SETTINGS_H

#include <vector>
#include <string>
#include "utils/floatpoint.h"
#include "utils\string.h"


//v1.04: modified timeEstimate for improving accuracy and added time information to GCode.
//v1.05: added arguements for home line
#define VERSION "V1.1.12"

#define FIX_HORRIBLE_UNION_ALL_TYPE_A    0x01
#define FIX_HORRIBLE_UNION_ALL_TYPE_B    0x02
#define FIX_HORRIBLE_EXTENSIVE_STITCHING 0x04
#define FIX_HORRIBLE_UNION_ALL_TYPE_C    0x08
#define FIX_HORRIBLE_KEEP_NONE_CLOSED    0x10
/**
 * Type of support material.
 * Grid is a X/Y grid with an outline, which is very strong, provides good support. But in some cases is hard to remove.
 * Lines give a row of lines which break off one at a time, making them easier to remove, but they do not support as good as the grid support.
 */
#define SUPPORT_TYPE_GRID                0
#define SUPPORT_TYPE_LINES               1
#define SUPPORT_TYPE_LINES_END_CONNECT	 2 //cubicon option 

#ifndef DEFAULT_CONFIG_PATH
#define DEFAULT_CONFIG_PATH "default.cfg"
#endif

#if defined(WIN32) || defined(WIN64)
  #define snprintf _snprintf
  #define vsnprintf _vsnprintf
  #define strcasecmp _stricmp
  #define strncasecmp _strnicmp
#else
  #define strcasecmp cura::stringcasecompare
#endif

/**
 * RepRap flavored GCode is Marlin/Sprinter/Repetier based GCode.
 *  This is the most commonly used GCode set.
 *  G0 for moves, G1 for extrusion.
 *  E values give mm of filament extrusion.
 *  Retraction is done on E values with G1. Start/end code is added.
 *  M106 Sxxx and M107 are used to turn the fan on/off.
 **/
#define GCODE_FLAVOR_REPRAP              0
/**
 * UltiGCode flavored is Marlin based GCode.
 *  UltiGCode uses less settings on the slicer and puts more settings in the firmware. This makes for more hardware/material independed GCode.
 *  G0 for moves, G1 for extrusion.
 *  E values give mm^3 of filament extrusion. Ignores the filament diameter setting.
 *  Retraction is done with G10 and G11. Retraction settings are ignored. G10 S1 is used for multi-extruder switch retraction.
 *  Start/end code is not added.
 *  M106 Sxxx and M107 are used to turn the fan on/off.
 **/
#define GCODE_FLAVOR_ULTIGCODE           1
/**
 * Makerbot flavored GCode.
 *  Looks a lot like RepRap GCode with a few changes. Requires MakerWare to convert to X3G files.
 *   Heating needs to be done with M104 Sxxx T0
 *   No G21 or G90
 *   Fan ON is M126 T0 (No fan strength control?)
 *   Fan OFF is M127 T0
 *   Homing is done with G162 X Y F2000
 **/
#define GCODE_FLAVOR_MAKERBOT           2

/**
cubicreator flavored is reprap based GCode.
**/
#define GCODE_FLAVOR_CUBICREATOR		3

#define MAX_EXTRUDERS 16

//#define DEBUG_FORCE_OPTIONS
//#define DEBUG_CURA_FORCE_OPTIONS

#define MINIMAL_POLYGON_LENGTH    10

#define RAFT_TYPE_CURA		0
#define RAFT_TYPE_CUBICON	1

class _ConfigSettingIndex
{
public:
    const char* key;
    int* ptr;

    _ConfigSettingIndex(const char* key, int* ptr) : key(key), ptr(ptr) {}
};

struct LayerControl
{
	int layerIndex = -1;
	std::string gcode;	

	LayerControl():gcode(){	}
	LayerControl(int _index, const char *pStr) :layerIndex(_index), gcode(pStr)
	{
	}
};
class ConfigSettings
{
private:
    std::vector<_ConfigSettingIndex> _index;
public:
    int layerThickness;
    int initialLayerThickness;
    int filamentDiameter;
    int filamentFlow;
    int extrusionWidth;
    int insetCount;
    int downSkinCount;
    int upSkinCount;
    int sparseInfillLineDistance;
    int infillOverlap;
    int skirtDistance;
    int skirtLineCount;
    int skirtMinLength;

    //Retraction settings
    int retractionAmount;
    int retractionAmountExtruderSwitch;
    int retractionSpeed;
    int retractionMinimalDistance;
    int minimalExtrusionBeforeRetraction;
    int retractionZHop;

	int retractRestoreCorrection;

    int enableCombing;
    int enableOozeShield;
    int wipeTowerSize;
    int multiVolumeOverlap;

    int initialSpeedupLayers;
    int initialLayerSpeed;
    int printSpeed;
    int infillSpeed;
    int inset0Speed;
    int insetXSpeed;
    int moveSpeed;
    int fanFullOnLayerNr;
	int skirtSpeed;
	int startSpeed;//cstyle

    //Support material
	int	supportInfillType;
    int supportAngle;
    int supportEverywhere;
    int supportLineDistance;
    int supportXYDistance;
    int supportZDistance;
    int supportExtruder;

    //Cool settings
    int minimalLayerTime;
    int minimalFeedrate;
    int coolHeadLift;
    int fanSpeedMin;
    int fanSpeedMax;

    //Raft settings
	int raftType;
    int raftMargin;
    int raftLineSpacing;
    int raftBaseThickness;
    int raftBaseLinewidth;
    int raftInterfaceThickness;
    int raftInterfaceLinewidth;

	int raftPlaneWidth1;
	int raftPlaneWidth2;
	int raftPlaneThickness1;
	int raftPlaneThickness2;
	int raftBaseSpeed;
	int raftInterfaceSpeed;
    
    FMatrix3x3 matrix;
    IntPoint objectPosition;
    int objectSink;

    int fixHorrible;
    int spiralizeMode;
    int gcodeFlavor;

	//additional option
	//cstyle add dash pattern
	int raftDashStride;
	int raftDashSize;
	int wallInOutOrder;
	int insetEndRetraction;
	int strongSupportLayerNr;
	int airGap;
	int layerSupport;
	int layerSupportTopMarginCount;
	int homelineThickness;

	//wipe head
	int			wipetype;//0, 1
	IntPoint	wipePosition;
	IntPoint	wipePosition1;
	int			wipeExtrusionAmount;

	int			finalizeRetraction;
	int			fanspeedThreshold;

    IntPoint extruderOffset[MAX_EXTRUDERS];
    const char* startCode;
    const char* endCode;

    ConfigSettings();
    bool setSetting(const char* key, const char* value);
    bool readSettings(void);
    bool readSettings(const char* path);
	std::vector<std::string> AnalyzeArgument(std::string str);

	//cstyle layercontrol test
	//bed temporature
	//ext temp
	//fan speed
	//pause
	std::vector<LayerControl> layerControls;

	//system setting
	int homeline_marginX;
	int homeline_marginY;
	int homeline_startX;
	int homeline_startY;
	int bed_width;
	int bed_height;
	
	int innerWallEndControlType;//0:endgap 1: endReduce
	int innerwallEndReduceRate;//0~100
	int innerWallEndGap;
	int outerWallEndGap;

	int pullEndOfWallTravelDistance;
	int pullEndOfWallAngle;//+-0~180
	int travelSpeedMatching;
	int outerWallStartOffset;
	int fixedWallBeginAngle;
	int finalizeExtLift;
	
	int initialLayerFillType;//0: line, 1:concentric

	//int disableForceLand;//0: land a model at zero  1:disable land a model
};

#endif//SETTINGS_H