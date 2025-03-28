#ifndef _MY_IO_H_
#define _MY_IO_H_

#include "config.h"

int readConfig(char* configName, Config* config);
// int readData(const char* fileName, int numX, int numY, double* dataArray);
int readData(double* field, int numX, int numY, const char* caseName, int timeStep, char* fieldName);
int writeData(double* field, int numX, int numY, const char* caseName, int timeStep, char* fieldName, int nx, int ny, int writePrec);
int writeDataHeader(char* fileName, double* field, int numX, int numY, char* header);

int write(char* filename, int dataNum, int num, char** headings, ...);

// directories
int createDirectoryIfNotExists(const char* dirPath);
int getNumericDirectories(const char* caseDirPath, double dirList[], int maxDirs, double* maxValue);

#endif