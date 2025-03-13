#include "stdio.h"

#include "myConst.h"
// #include "config.h"
#include "myIO.h"

int main()
{
    const char* caseName = "test_01";
    
    Config cfg;
    char cfgName[128] = "";
    strcpy(cfgName,caseName); // cfgName = caseName
    strcat(cfgName,"/config"); // cfgName = caseName/config

    readConfig(cfgName,&cfg);

    const int nx = cfg.nx;
    const int ny = cfg.ny;

    double* u = (double*)calloc((nx+3)*(ny+2), sizeof(double));

    readData(u,nx+3,ny+2,caseName,0,"U");

    printf("%lf\n", u[1 + 3*(nx+3)]);


    return 0;
}