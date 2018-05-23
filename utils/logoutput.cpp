/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <stdio.h>
#include <stdarg.h>

#include "logoutput.h"

static int verbose_level;
static bool progressLogging;
static ANSWERCB pCallBack = NULL;
#define PROGRESS_FMT "Progress:%s:%i:%i\n"
#define LOG_BUFF_SIZE (4096)

void setLogCallBack(ANSWERCB cb)
{
	pCallBack = cb;
}

void increaseVerboseLevel()
{
    verbose_level++;
}

void enableProgressLogging()
{
    progressLogging = true;
}

void logError(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
	if (NULL != pCallBack)
	{
		char buff[LOG_BUFF_SIZE] = { 0, };
		vsprintf(buff, fmt, args);
		//vfprintf(buff, fmt, args);
		pCallBack(buff);
	}
	va_end(args);
	fflush(stdout);
}

void log(const char* fmt, ...)
{
    if (verbose_level < 1)
        return;

    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
	if (NULL != pCallBack)
	{
		char buff[LOG_BUFF_SIZE] = { 0, };
		vsprintf(buff, fmt, args);
		pCallBack(buff);
	}
	va_end(args);
	fflush(stdout);
}

void logProgress(const char* type, int value, int maxValue)
{
    if (!progressLogging)
        return;

	if (NULL != pCallBack)
	{
		char buff[LOG_BUFF_SIZE] = { 0, };
		sprintf(buff, PROGRESS_FMT, type, value, maxValue);
		pCallBack(buff);
	}
	else
	{
		fprintf(stdout, PROGRESS_FMT, type, value, maxValue);
		fflush(stdout);
	}
}
