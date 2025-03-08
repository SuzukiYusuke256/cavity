#include <stdio.h>
#include <stdarg.h> 
#include <stdlib.h>

int writeData(char* filename, int dataNum, int num, char** headings, ...);

int main()
{
    char filename[100] = "output";
    int num = 5;
    int dataNum = 3;
    double** dataPointerList = (double**)calloc(dataNum, sizeof(double*));

    double* data1 = (double*)calloc(num, sizeof(double));
    double* data2 = (double*)calloc(num, sizeof(double));
    double* data3 = (double*)calloc(num, sizeof(double));

    // save the pointer of data in one array 
    // so that it will be easier to free memory at the end of the code
    dataPointerList[0] = data1;
    dataPointerList[1] = data2;
    dataPointerList[2] = data3;

    data1[2] = 256;

    char *dataNameArray[] = {
        "Apple",
        "Banana",
        "Cherry"
    };

    // write data
    // write(filename, 3, num, dataNameArray, data1, data2, data3);
    writeData(filename, dataNum, num, dataNameArray, data1, data2, data3);

    // free memory
    for(int i=0; i<dataNum; i++)
    {
        free(dataPointerList[i]);
    }
    free(dataPointerList);
}

// write data
// data must be array of double
int writeData(char* filename, int dataNum, int num, char** headings, ...)
{
    // open file
    FILE *file = fopen(filename, "w");

    // file open error check
    if (file == NULL) 
    {
        printf("Failed to open file\n");
        return 1;
    }
    
    // variable array
    va_list args;
    va_start(args, headings);

    // write headings
    fprintf(file, "Number");
    for (int i = 0; i < dataNum; i++)
    {
        fprintf(file, " %s",headings[i]);
    }
    fprintf(file, "\n");

    // obtain the pointer of each data
    double **dataPointerList = (double**)calloc(dataNum, sizeof(double*));

    for (int i = 0; i < dataNum; i++)
    {
        dataPointerList[i] = va_arg(args, double*);  // obtain next pointer
    }

    // write data
    for (int i = 0; i < num; i++)
    {
        // write time step number
        fprintf(file, "%d", i);

        // write each data 
        for (int j = 0; j < dataNum; j++)
        {
            fprintf(file, " %e", dataPointerList[j][i]);
        }

        // add new line
        fprintf(file, "\n");
    }

    va_end(args);

    // close file
    fclose(file);

    // free memory
    free(dataPointerList);

    return 0;
}

