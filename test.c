#include "stdio.h"

#include "myConst.h"
// #include "config.h"
#include "myIO.h"

int main()
{
    // double* configArray = (double*)calloc(CFG_NUM,sizeof(double));
    Config config;
    readConfig("test_02/config",&config);

    printf("%s\n",config.caseName);
    printf("%.3e\n",config.convergenceThreshold);

    return 0;
}