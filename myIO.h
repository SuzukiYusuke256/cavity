#ifndef _MY_IO_H_
#define _MY_IO_H_

int readData(const char* fileName, int numX, int numY, double* dataArray);
int writeData(char* fileName, double* field, int numX, int numY, char* header);
int write(char* filename, int dataNum, int num, char** headings, ...);

#endif