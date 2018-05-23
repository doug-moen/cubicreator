/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef EXPORT_DLL_LIBRARY
#include <windows.h>
#define DLL_LIBRARY (1)
#else
#include <sys/time.h>
#define DLL_LIBRARY (0)
#endif

#include <signal.h>
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#include <execinfo.h>
#include <sys/resource.h>
#endif
#include <stddef.h>

#include "utils/gettime.h"
#include "utils/logoutput.h"
#include "sliceDataStorage.h"

#include "modelFile/modelFile.h"
#include "settings.h"
#include "optimizedModel.h"
#include "multiVolumes.h"
#include "polygonOptimizer.h"
#include "slicer.h"
#include "layerPart.h"
#include "inset.h"
#include "skin.h"
#include "infill.h"
#include "bridge.h"
#include "support.h"
#include "pathOrderOptimizer.h"
#include "skirt.h"
#include "raft.h"
#include "comb.h"
#include "gcodeExport.h"
#include "fffProcessor.h"
#include "utils/util.h"

void print_usage()
{
    //log("usage: CuraEngine [-h] [-v] [-m 3x3matrix] [-c <config file>] [-s <settingkey>=<value>] -o <output.gcode> <model.stl>\n");
	log("usage: CubiEngine [-h] [-v] [-m 3x3matrix] [-c <config file>] [-s <settingkey>=<value>] -o <output.gcode> <model.stl>\n");
}


#if DLL_LIBRARY
SliceDataStorage *g_SliceData = NULL;
#if WIN64
__declspec(dllexport) int SliceModel64(SModel *pModel, float sliceZ, bool dirty, bool keepNoneClosed, bool extensiveStitching)
#else
__declspec(dllexport) int SliceModel32(SModel *pModel, float sliceZ, bool dirty, bool keepNoneClosed, bool extensiveStitching)
#endif
{
	if (pModel)
	{
		pModel->reserved += 1;
		//for (int i = 0; i < pModel->size; i += 3)
		//{
		//	float fX = pModel->pBuff[i + 0];
		//	float fY = pModel->pBuff[i + 1];
		//	float fZ = pModel->pBuff[i + 2];
		//	logError("[%d]:[%f][%f][%f]", i / 3, fX, fY, fZ);
		//}
		if (dirty && NULL != g_SliceData)
		{
			delete g_SliceData;
			g_SliceData = NULL;
		}

		g_SliceData = Slicer::Slice(g_SliceData, pModel, sliceZ, keepNoneClosed, extensiveStitching);
		return 0;
	}
	return -1;
}
#endif

//Signal handler for a "floating point exception", which can also be integer division by zero errors.
void signal_FPE(int n)
{
    (void)n;
    logError("Arithmetic exception.\n");
    exit(1);
}

#if DLL_LIBRARY && WIN64
__declspec(dllexport) int CubiEngine64(int argc, char **argv, ANSWERCB cb)
#elif DLL_LIBRARY && WIN32
__declspec(dllexport) int CubiEngine32(int argc, char **argv, ANSWERCB cb)
#else
int CubiEngine(int argc, char **argv, ANSWERCB cb)
#endif
{
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
	//Lower the process priority on linux and mac. On windows this is done on process creation from the GUI.
	setpriority(PRIO_PROCESS, 0, 10);
#endif

	setLogCallBack(cb);
	//Register the exception handling for arithmic exceptions, this prevents the "something went wrong" dialog on windows to pop up on a division by zero.
	signal(SIGFPE, signal_FPE);

	ConfigSettings config;
	fffProcessor processor(config);

	logError("CubiEngine version %s\n", VERSION);

	if (!config.readSettings()) {
		logError("Default config '%s' not used\n", DEFAULT_CONFIG_PATH);
	}

	enableProgressLogging();

	//std::vector<char*> test;
	//test.push_back(argv[argn]);

#if DLL_LIBRARY
	for (int argn = 0; argn < argc; argn++)
#else
	for (int argn = 1; argn < argc; argn++)
#endif
	{
		char* str = argv[argn];
		if (str[0] == '-')
		{
			for (str++; *str; str++)
			{
				switch (*str)
				{
				case 'h':
					print_usage();
					exit(1);
				case 'v':
					increaseVerboseLevel();
					break;
				case 'p':
					enableProgressLogging();
					break;
				case 'g':
					argn++;
					//Connect the GUI socket to the given port number.
					processor.guiConnect(atoi(argv[argn]));
					break;
				case 'b':
					argn++;
					//The binaryMeshBlob is depricated and will be removed in the future.
					binaryMeshBlob = fopen(argv[argn], "rb");
					break;
				case 'o':
					argn++;
					if (!processor.setTargetFile(argv[argn]))
					{
						logError("Failed to open %s for output.\n", argv[argn]);
						exit(1);
					}
					break;
				case 'c':
				{
							// Read a config file from the given path
							argn++;
							if (!config.readSettings(argv[argn])) {
								logError("Failed to read config '%s'\n", argv[argn]);
							}
				}
					break;
				case 's':
				{
							//Parse the given setting and store it.
							argn++;
							char* valuePtr = strchr(argv[argn], '=');
							if (valuePtr)
							{
								*valuePtr++ = '\0';

								if (!config.setSetting(argv[argn], valuePtr))
									logError("Setting not found: %s %s\n", argv[argn], valuePtr);
							}
				}
					break;
				case 'm':
					//Read the given rotation/scale matrix
					argn++;
					sscanf(argv[argn], "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
						&config.matrix.m[0][0], &config.matrix.m[0][1], &config.matrix.m[0][2],
						&config.matrix.m[1][0], &config.matrix.m[1][1], &config.matrix.m[1][2],
						&config.matrix.m[2][0], &config.matrix.m[2][1], &config.matrix.m[2][2]);
					break;
				default:
					logError("Unknown option: %c\n", *str);
					break;
				}
			}
		}
		else{
			try {
				//Catch all exceptions, this prevents the "something went wrong" dialog on windows to pop up on a thrown exception.
				// Only ClipperLib currently throws exceptions. And only in case that it makes an internal error.
				processor.processFile(argv[argn]);
			}
			catch (...){
				logError("Unknown exception\n");
				exit(1);
			}
		}
	}

	//Finalize the processor, this adds the end.gcode. And reports statistics.
	processor.finalize();	

	if (NULL != binaryMeshBlob)
		fclose(binaryMeshBlob);
	return 0;
}

#if !DLL_LIBRARY
int main(int argc, char **argv)
{
	return CubiEngine(argc, argv, NULL);
}
#endif
