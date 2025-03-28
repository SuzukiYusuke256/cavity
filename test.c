#include "stdio.h"

#include "myConst.h"
#include "config.h"
#include "myIO.h"

int main()
{
    const char* caseDirPath = "cavity_smac_01";
    double dirList[1000] = {0};
    int maxDirs = 1000;
    double maxValue = 0.0;

    const int count = getNumericDirectories(caseDirPath, dirList, maxDirs, &maxValue);

    printf("Number of numeric directories: %d\n", count);
    printf("Max number: %lf\n", maxValue);

    return 0;
}