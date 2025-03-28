// cavity flow solver using SMAC method
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "myIO.h"
#include "myConst.h"
#include "config.h"

// define functions
void correctBoundaryConditions(double* u, double* v, double* p, const int* uBC, const int* vBC, const int* pBC, 
                               const double* uBCVal, const double* vBCVal, const double* pBCVal, int xn, int yn, double dx, double dy);
void predictVelocity(double* u, double* v, double* p, double* du, double* dv, double* dp, double Re, double dx, double dy, double dt, int xn, int yn);
void poissonPressure(double* u, double* v, double* p, double* du, double* dv, double* dp, double Re, double dx, double dy, double dt, int xn, int yn);
void correct(double* u, double* v, double* p, double* du, double* dv, double* dp, double dx, double dy, double dt, int xn, int yn, int pRefID);
void writeAll(int timeStep, double* u, double* v, double* p, double* du, double* dv, double* dp, int xn, int yn, char* dataDirName);
void writeLog(char* fileName, int timeStep, double dt, time_t elapsedTime, double* res, int resNum);
void calcResidual(double* field, int startIndexX, int startIndexY, int numX, int numY, int maxNumX, double* res, int resStartIndex);

int main(int argc, char* argv[])
{
    // コマンドライン引数のチェック
    if (argc != 2) {
        printf("Usage: %s <caseName>\n", argv[0]);
        return 1;
    }

    // コマンドライン引数からcaseNameを取得
    const char* caseName = argv[1];
    // const char* caseName = "test_03";
    
    Config cfg;
    char cfgName[128] = "";
    strcpy(cfgName,caseName); // cfgName = caseName
    strcat(cfgName,"/config"); // cfgName = caseName/config
    readConfig(cfgName,&cfg);

    // if the case name does not match the config file, return error
    if (strcmp(caseName, cfg.caseName) != 0) 
    {
        fprintf(stderr, "Error: caseName does not match cfg.caseName\n");
        return -1;
    }

    const int stepNum           = cfg.stepNum;
    const int outputInterval    = cfg.outputInterval;
    const double dt             = cfg.deltaT;
    const int xn                = cfg.nx;
    const int yn                = cfg.ny;
    const double Re             = cfg.re;
    const int withIC            = cfg.withInitialCondition;
    const double thresh         = cfg.convergenceThreshold;
    const int writePrec         = cfg.writePrecision; // precision for output data

    const int logOutputInterval = outputInterval/10; // ログの出力間隔

    // 圧力の基準セル．この格子点の圧力を0とする．
    const int pRefID = 1*(xn+2) + 1; // i = 1, j = 1 (0始まり)  

    // boundary conditions /////////////////////////////////////////////////

    // couette flow 
    // const int uBC[4] = {0,0,1,1}; // uの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    // const int vBC[4] = {0,0,1,1}; // vの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    // const int pBC[4] = {1,1,1,1}; // pの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann

    // const double uBCVal[4] = {1.0, 0.0, 0.0, 0.0}; // uの境界条件値 上下左右の順
    // const double vBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // vの境界条件値 上下左右の順
    // const double pBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // pの境界条件値 上下左右の順

    // 2D poiseuille flow
    // const int uBC[4] = {0,0,1,1}; // uの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    // const int vBC[4] = {0,0,1,1}; // vの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    // const int pBC[4] = {1,1,1,1}; // pの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann

    // const double uBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // uの境界条件値 上下左右の順
    // const double vBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // vの境界条件値 上下左右の順
    // const double pBCVal[4] = {0.0, 0.0, -1.0, -1.0}; // pの境界条件値 上下左右の順

    // cavity flow 
    const int uBC[4] = {0,0,0,0}; // uの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    const int vBC[4] = {0,0,0,0}; // vの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann
    const int pBC[4] = {1,1,1,1}; // pの境界条件タイプ 上下左右の順 0: Dirichlet 1: Neumann

    const double uBCVal[4] = {1.0, 0.0, 0.0, 0.0}; // uの境界条件値 上下左右の順
    const double vBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // vの境界条件値 上下左右の順
    const double pBCVal[4] = {0.0, 0.0, 0.0, 0.0}; // pの境界条件値 上下左右の順
    
    // output file name
    char logFileName[100] = "log"; // default 
    char timeDirName[100] = {0}; 
    char tmpDirName[100] = {0}; 
    char header[1000] = {0}; // header for data output

    sprintf(logFileName,"%s/log",caseName);

    // 計算時間 計測
    const time_t startTime = time(NULL); // プログラムの開始時間を取得
    time_t elapsedTime = time(NULL) - startTime; // 経過時間．ログに出力

    double dx = 1.0 / (double)xn;
    double dy = 1.0 / (double)yn;
    
    // field initialization
    double* u  = (double*)calloc((xn+3)*(yn+2),sizeof(double));
    double* v  = (double*)calloc((xn+2)*(yn+3),sizeof(double));
    double* p  = (double*)calloc((xn+2)*(yn+2),sizeof(double));
    double* du = (double*)calloc((xn+3)*(yn+2),sizeof(double));
    double* dv = (double*)calloc((xn+2)*(yn+3),sizeof(double));
    double* dp = (double*)calloc((xn+2)*(yn+2),sizeof(double));

    double residuals[6] = {0}; // ログに残差を出力するための配列．0で初期化

    // initial conditions
    if (withIC == 1)
    {   
        printf("reading initial conditions...\n");
        // readData("initialConditions/U",xn+3,yn+2,u);
        // readData("initialConditions/V",xn+2,yn+3,v);
        // readData("initialConditions/p",xn+2,yn+2,p);
    }
    else{

        // boundary conditions
        correctBoundaryConditions(u,v,p,uBC,vBC,pBC,uBCVal,vBCVal,pBCVal,xn,yn,dx,dy);

        // output directory for computed data
        sprintf(timeDirName,"%s/%d",caseName,0);
        createDirectoryIfNotExists(timeDirName);

        // write initial conditions
        writeData(u,xn+3,yn+2,caseName,0,"U",xn,yn,writePrec);
        writeData(v,xn+2,yn+3,caseName,0,"V",xn,yn,writePrec);
        writeData(p,xn+2,yn+2,caseName,0,"p",xn,yn,writePrec);
    }

    // ログファイルのヘッダを書き込み
    FILE* outputLog = fopen(logFileName,"w");

    // ファイルが正しく開けたかを確認
    if (outputLog == NULL) {
        printf("Failed to open file\n");
        return 0;
    }
    else
    {
        fprintf(outputLog,"TimeStep Time URes URes(cell) VRes VRes(cell) pRes pRes(cell) executionTime[s]\n"); // ヘッダを書き込み
        fclose(outputLog); // ファイルを閉じる
    }

    // SMAC法での計算
    for(int kk=0; kk<=stepNum; kk++)
    {

        // 仮速度の計算
        // u,vに仮速度u*,v*が入る
        predictVelocity(u,v,p,du,dv,dp,Re,dx,dy,dt,xn,yn);

        // 圧力ポアソン方程式を解く
        // 圧力修正子を計算
        poissonPressure(u,v,p,du,dv,dp,Re,dx,dy,dt,xn,yn);

        // 速度，圧力の修正
        correct(u,v,p,du,dv,dp,dx,dy,dt,xn,yn,pRefID);

        // 境界条件の修正
        correctBoundaryConditions(u,v,p,uBC,vBC,pBC,uBCVal,vBCVal,pBCVal,xn,yn,dx,dy);

        // 残差の計算
        calcResidual(du,2,1,xn-1,yn,xn+3,residuals,0);
        calcResidual(dv,1,2,xn,yn-1,xn+2,residuals,2);
        calcResidual(dp,1,1,xn,yn,xn+2,residuals,4);

        // 計算の終了判定
        int isEnd = 0;

        // 最大回数に達する
        if (kk == stepNum)
        {
            isEnd = 1;
        }

        // 残差が収束判定のしきい値を下回る kk=0は残差が0なので除く
        if (kk>0 && residuals[1] < thresh && residuals[3] < thresh) // du, dvで収束を判定
        {
            isEnd = 1; 
        } 

        // kk=0のデータは初期条件なので出力しない．
        if (kk == 0)
        {
            continue;
        }

        // 残差をログに出力．計算終了時は別でログを出力するため，終了条件を満たしていない場合のみ出力
        if (kk % logOutputInterval == 0 && isEnd == 0)
        {
            // ログに残差を出力
            elapsedTime = time(NULL) - startTime;
            writeLog(logFileName,kk,dt,elapsedTime,residuals,6); 
        }

        // 残差を出力. 一定間隔かつ収束条件を満たしていない場合．
        if (kk % outputInterval == 0 && isEnd == 0)
        {
            // output directory for computed data
            sprintf(timeDirName,"%s/%d",caseName,kk);
            createDirectoryIfNotExists(timeDirName);

            writeData(u,xn+3,yn+2,caseName,kk,"U",xn,yn,writePrec);
            writeData(du,xn+3,yn+2,caseName,kk,"dU",xn,yn,writePrec);
            writeData(v,xn+2,yn+3,caseName,kk,"V",xn,yn,writePrec);
            writeData(dv,xn+2,yn+3,caseName,kk,"dV",xn,yn,writePrec);
            writeData(p,xn+2,yn+2,caseName,kk,"p",xn,yn,writePrec);
            writeData(dp,xn+2,yn+2,caseName,kk,"dp",xn,yn,writePrec);

        }

        // 収束条件を満たした場合，計算を終了
        if (isEnd == 1)
        {
            // ログに残差を出力
            elapsedTime = time(NULL) - startTime;
            writeLog(logFileName,kk,dt,elapsedTime,residuals,6); 
            
            // output directory for computed data
            sprintf(timeDirName,"%s/%d",caseName,kk);
            createDirectoryIfNotExists(timeDirName);

            writeData(u,xn+3,yn+2,caseName,kk,"U",xn,yn,writePrec);
            writeData(du,xn+3,yn+2,caseName,kk,"dU",xn,yn,writePrec);
            writeData(v,xn+2,yn+3,caseName,kk,"V",xn,yn,writePrec);
            writeData(dv,xn+2,yn+3,caseName,kk,"dV",xn,yn,writePrec);
            writeData(p,xn+2,yn+2,caseName,kk,"p",xn,yn,writePrec);
            writeData(dp,xn+2,yn+2,caseName,kk,"dp",xn,yn,writePrec);

            // 計算終了を表示
            printf("%s finished\n",caseName);

            break;
        }
    }

    // メモリの解放
    free(u);
    free(v);
    free(p);
    free(du);
    free(dv);
    free(dp);

    // printf("%lf\n",u[3]);

}

// 境界条件に合わせてセルの値を更新
void correctBoundaryConditions(double* u, double* v, double* p, const int* uBC, const int* vBC, const int* pBC, 
    const double* uBCVal, const double* vBCVal, const double* pBCVal, int xn, int yn, double dx, double dy)
{
    // U boundary conditions ////////////////////////////////////////////////////////////

    // top
    if (uBC[0] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+2; i++)
        {
            u[i + (yn+1) * (xn+3)] = 2.0 * uBCVal[0] - u[i + yn * (xn+3)];
        }
    }
    else if (uBC[0] == 1) // neumann condition
    {
        for(int i=1; i<xn+2; i++)
        {
            u[i + (yn+1) * (xn+3)] = u[i + yn * (xn+3)] + uBCVal[0] * dy; // gradient perpendicular to the wall is specified.
        }
    }

    // bottom
    if (uBC[1] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+2; i++)
        {
            u[i + 0 * (xn+3)] = 2.0 * uBCVal[1] - u[i + 1 * (xn+3)];
        }
    }
    else if (uBC[1] == 1) // neumann condition
    {
        for(int i=1; i<xn+2; i++)
        {
            u[i + 0 * (xn+3)] = u[i + 1 * (xn+3)] - uBCVal[1] * dy; // gradient perpendicular to the wall is specified.
        }
    }

    // left
    if (uBC[2] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+1; j++)
        {
            u[1 + j * (xn+3)] = uBCVal[2]; // boudary surface velocity
            u[0 + j * (xn+3)] = u[2 + j * (xn+3)]; // boundary condition for continuity
        }
    }
    else if (uBC[2] == 1) // neumann condition
    {
        for(int j=1; j<yn+1; j++)
        {
            u[1 + j * (xn+3)] = u[2 + j * (xn+3)] + uBCVal[2] * dx; // gradient is constant perpendiculer to the wall
            u[0 + j * (xn+3)] = u[2 + j * (xn+3)] + uBCVal[2] * 2.0 * dx; 
        }
    }

    // right
    if (uBC[3] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+1; j++)
        {
            u[(xn+1) + j * (xn+3)] = uBCVal[3]; // boudary surface velocity
            u[(xn+2) + j * (xn+3)] = u[xn + j * (xn+3)]; // boundary condition for continuity
        }
    }
    else if (uBC[3] == 1) // neumann condition
    {
        for(int j=1; j<yn+1; j++)
        {
            u[(xn+1) + j * (xn+3)] = u[xn + j * (xn+3)] + uBCVal[3] * dx; 
            u[(xn+2) + j * (xn+3)] = u[xn + j * (xn+3)] + uBCVal[3] * 2.0 * dx;
        }
    }

    // V boundary conditions ////////////////////////////////////////////////////////////

    // top
    if (vBC[0] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+1; i++)
        {
            v[i + (yn+1) * (xn+2)] = vBCVal[0]; // boudary surface velocity
            v[i + (yn+2) * (xn+2)] = v[i + yn * (xn+2)]; // boundary condition for continuity
        }
    }
    else if (vBC[0] == 1) // neumann condition
    {
        for(int i=1; i<xn+1; i++)
        {
            v[i + (yn+1) * (xn+2)] = v[i + yn * (xn+2)] + vBCVal[0] * dy; // gradient is constant perpendiculer to the wall
            v[i + (yn+2) * (xn+2)] = v[i + yn * (xn+2)] + vBCVal[0] * 2.0 * dy; 
        }
    }

    // bottom
    if (vBC[1] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+1; i++)
        {
            v[i + 1 * (xn+2)] = vBCVal[1]; // boudary surface velocity
            v[i + 0 * (xn+2)] = v[i + 2 * (xn+2)]; // boundary condition for continuity
            
            // u[1 + j * (xn+3)] = uBCVal[2]; // boudary surface velocity
            // u[0 + j * (xn+3)] = u[2 + j * (xn+3)]; // boundary condition for continuity
        }
    }
    else if (vBC[1] == 1) // neumann condition
    {
        for(int i=1; i<xn+1; i++)
        {
            v[i + 1 * (xn+2)] = v[i + 2 * (xn+2)] + vBCVal[1] * dy; // gradient is constant perpendiculer to the wall
            v[i + 0 * (xn+2)] = v[i + 2 * (xn+2)] + vBCVal[1] * 2.0 * dy; 
        }
    }

    // left
    if (vBC[2] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+2; j++)
        {
            v[0 + j * (xn+2)] = 2.0 * vBCVal[2] - v[1 + j * (xn+2)];
        }
    }
    else if (vBC[2] == 1) // neumann condition
    {
        for(int j=1; j<yn+2; j++)
        {
            v[0 + j * (xn+2)] = v[1 + j * (xn+2)] - vBCVal[2] * dx; // gradient perpendicular to the wall is specified.
        }
    }

    // right
    if (vBC[3] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+2; j++)
        {
            v[(xn+1) + j * (xn+2)] = 2.0 * vBCVal[3] - v[xn + j * (xn+2)];
        }
    }
    else if (vBC[3] == 1) // neumann condition
    {
        for(int j=1; j<yn+2; j++)
        {
            v[(xn+1) + j * (xn+2)] = v[xn + j * (xn+2)] + vBCVal[3] * dx; // gradient perpendicular to the wall is specified.
        }
    }

    // pressure boundary conditions
    // top
    if (pBC[0] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+1; i++)
        {
            p[i + (yn+1) * (xn+2)] = 2.0 * pBCVal[0] - p[i + yn * (xn+2)];
        }
    }
    else if (pBC[0] == 1) // neumann condition
    {
        for(int i=1; i<xn+1; i++)
        {
            p[i + (yn+1) * (xn+2)] = p[i + yn * (xn+2)] + pBCVal[0] * dy; // gradient perpendicular to the wall is specified.
        }
    }

    // bottom
    if (pBC[1] == 0) // dirichlet condition
    {
        for(int i=1; i<xn+1; i++)
        {
            p[i + 0 * (xn+2)] = 2.0 * pBCVal[1] - p[i + 1 * (xn+2)];
        }
    }
    else if (pBC[1] == 1) // neumann condition
    {
        for(int i=1; i<xn+1; i++)
        {
            p[i + 0 * (xn+2)] = p[i + 1 * (xn+2)] - pBCVal[1] * dy; // gradient perpendicular to the wall is specified.
        }
    }

    // left
    if (pBC[2] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+1; j++)
        {
            p[0 + j * (xn+2)] = 2.0 * pBCVal[2] - p[1 + j * (xn+2)];
        }
    }
    else if (pBC[2] == 1) // neumann condition
    {
        for(int j=1; j<yn+1; j++)
        {
            p[0 + j * (xn+2)] = p[1 + j * (xn+2)] - pBCVal[2] * dx; // gradient perpendicular to the wall is specified.
        }
    }

    // right
    if (pBC[3] == 0) // dirichlet condition
    {
        for(int j=1; j<yn+1; j++)
        {
            p[(xn+1) + j * (xn+2)] = 2.0 * pBCVal[3] - p[xn + j * (xn+2)];
        }
    }
    else if (pBC[3] == 1) // neumann condition
    {
        for(int j=1; j<yn+1; j++)
        {
            p[(xn+1) + j * (xn+2)] = p[xn + j * (xn+2)] + pBCVal[3] * dx; // gradient perpendicular to the wall is specified.
        }
    }
}

void predictVelocity(double* u, double* v, double* p, double* du, double* dv, double* dp, double Re, double dx, double dy, double dt, int xn, int yn)
{
    // NS方程式で仮速度を計算する
    double ududx = 0;
    double vdudy = 0;
    double udvdx = 0;
    double vdvdy = 0;

    double dpdx = 0;
    double dpdy = 0;

    double ddudxx = 0;
    double ddudyy = 0;
    double ddvdxx = 0;
    double ddvdyy = 0;

    // U
    for(int j=1; j<yn+1; j++)
    {
        for(int i=1; i<xn+2; i++)
        {
            ududx = u[i + j*(xn+3)] * (u[(i+1) + j *(xn+3)] - u[(i-1) + j * (xn+3)]) / (2.0 * dx);
            vdudy = (v[(i-1) + j * (xn+2)] + v[i + j * (xn+2)] + v[(i-1) + (j+1) * (xn+2)] + v[i + (j+1) * (xn+2)]) / 4.0 * (u[i +(j+1)*(xn+3)] - u[i+(j-1)*(xn+3)]) / (2.0 * dy);

            dpdx = (p[i + j*(xn+2)] - p[(i-1) + j*(xn+2)]) / dx;

            ddudxx = (u[(i+1) + j*(xn+3)] - 2*u[i + j*(xn+3)] + u[(i-1) + j*(xn+3)]) / (dx * dx);
            ddudyy = (u[i + (j+1)*(xn+3)] - 2*u[i + j*(xn+3)] + u[i + (j-1)*(xn+3)]) / (dy * dy);

            du[i + j*(xn+3)] = dt * (-(ududx + vdudy) - dpdx + (ddudxx + ddudyy) / Re);
            // u[j + i*(xn+3)] += dt * (-(ududx + vdudy) - dpdx + (ddudxx + ddudyy) / Re);
        }
    }

    // V
    for(int j=1; j<yn+2; j++)
    {
        for(int i=1; i<xn+1; i++)
        {
            udvdx = (u[(i+1) + j * (xn+3)] + u[i + j * (xn+3)] + u[(i+1) + (j-1) * (xn+3)] + u[i + (j-1) * (xn+3)]) / 4.0 * (v[(i+1) + j * (xn+2)] - v[(i-1) + j * (xn+2)]) / (2.0 * dx);
            vdvdy = v[i + j * (xn+2)] * (v[i + (j+1) * (xn+2)] - v[i + (j-1) * (xn+2)]) / (2.0 * dy);

            dpdy = (p[i + j*(xn+2)] - p[i + (j-1)*(xn+2)]) / dy;

            ddvdxx = (v[(i+1) + j*(xn+2)] - 2*v[i + j*(xn+2)] + v[(i-1) + j*(xn+2)]) / (dx * dx);
            ddvdyy = (v[i + (j+1)*(xn+2)] - 2*v[i + j*(xn+2)] + v[i + (j-1)*(xn+2)]) / (dy * dy);

            dv[i + j*(xn+2)] = dt * (-(udvdx + vdvdy) - dpdy + (ddvdxx + ddvdyy) / Re);
            // v[j + i*(xn+2)] += dt * (-(udvdx + vdvdy) - dpdy + (ddvdxx + ddvdyy) / Re);
        }
    }

    // 速度場を更新
    // U
    for(int j=1; j<yn+1; j++)
    {
        for(int i=1; i<xn+2; i++)
        {
            u[i + j*(xn+3)] = u[i + j*(xn+3)] + du[i + j*(xn+3)];
        }
    }
    // V
    for(int j=1; j<yn+2; j++)
    {
        for(int i=1; i<xn+1; i++)
        {
            v[i + j*(xn+2)] = v[i + j*(xn+2)] + dv[i + j*(xn+2)];
        }
    }

}

void poissonPressure(double* u, double* v, double* p, double* du, double* dv, double* dp, double Re, double dx, double dy, double dt, int xn, int yn)
{

    // 圧力修正量の計算
    double prevdp;
    double err;
    int itrMax = 1000;

    // dpを初期化
    for(int i=0; i<yn+2; i++)
    {
        for(int j=0; j<xn+2; j++)
        {
            dp[i+j*(xn+2)] = 0;
        }
    }

    for (int itr=0; itr<itrMax; itr++)
    {
        // 圧力境界条件
        // 壁面での勾配0
        // 上下面
        for(int i=0; i<xn+2; i++)
        {
            dp[i + 0*(xn+2)] = dp[i + 1*(xn+2)];
            dp[i + (yn+1)*(xn+2)] = dp[i + yn*(xn+2)];
        }
        //左右面
        for(int j=0; j<yn+2; j++)
        {
            dp[0 + j * (xn+2)] = dp[1 + j * (xn+2)];
            dp[(xn+1) + j * (xn+2)] = dp[xn + j * (xn+2)];
        }

        err = 0.0;
        
        // poisson方程式を解く
        for(int j=1; j<yn+1; j++)
        {
            for(int i=1; i<xn+1; i++)
            {
                prevdp = dp[i+j*(xn+2)];
                dp[i+j*(xn+2)] = 0.5 * (dx*dx*dy*dy) / (dx*dx + dy*dy) * ((dp[(i+1)+j*(xn+2)] + dp[(i-1)+j*(xn+2)]) / (dx*dx) + (dp[i+(j+1)*(xn+2)] + dp[i+(j-1)*(xn+2)]) / (dy*dy) - 1.0 / dt * ((u[(i+1)+j*(xn+3)] - u[i+j*(xn+3)]) / dx + (v[i+(j+1)*(xn+2)] - v[i+j*(xn+2)]) / dy));
                err = err + fabs(dp[i+j*(xn+2)] - prevdp);
                // printf("dp : %lf, prev : %lf\n",dp[j+i*(xn+2)],prevdp);
            }
        }

        if(err < 0.000001)
        {
            // printf("itr : %d, err : %lf\n",itr,err);
            break;
        }
    }
}

// 速度，圧力の修正
void correct(double* u, double* v, double* p, double* du, double* dv, double* dp, double dx, double dy, double dt, int xn, int yn, int pRefID)
{
    // 圧力と速度の修正
    // U
    for(int i=1; i<yn+1; i++)
    {
        for(int j=1; j<xn+2; j++)
        { 
            u[j+i*(xn+3)] += (-dt) * (dp[j+i*(xn+2)] - dp[(j-1)+i*(xn+2)]) / dx;
            du[j+i*(xn+3)] += (-dt) * (dp[j+i*(xn+2)] - dp[(j-1)+i*(xn+2)]) / dx; // 残差の計算用
        }
    }

    // V
    for(int i=1; i<yn+2; i++)
    {
        for(int j=1; j<xn+1; j++)
        {
            v[j+i*(xn+2)] += (-dt) * (dp[j+i*(xn+2)] - dp[j+(i-1)*(xn+2)]) / dy;
            dv[j+i*(xn+2)] += (-dt) * (dp[j+i*(xn+2)] - dp[j+(i-1)*(xn+2)]) / dy;
        }
    }


    // 圧力
    // 圧力修正子により圧力を修正．
    for(int i=1; i<yn+1; i++)
    {
        for(int j=1; j<xn+1; j++)
        {
            p[j+i*(xn+2)] += dp[j+i*(xn+2)];
        }
    }
    // 基準セルの圧力が0になるように全体を修正．
    double pRef = p[pRefID];
    for(int i=1; i<yn+1; i++)
    {
        for(int j=1; j<xn+1; j++)
        {
            p[j+i*(xn+2)] -= pRef;
        }
    }


}

// データを書き込む
void writeAll(int timeStep, double* u, double* v, double* p, double* du, double* dv, double* dp, int xn, int yn, char* dataDirName)
{

    // 出力
    FILE* Ures = fopen("data/U","w");
    FILE* Vres = fopen("data/V","w");
    FILE* Pres = fopen("data/p","w");
    FILE* dUres = fopen("data/dU","w");
    FILE* dVres = fopen("data/dV","w");
    FILE* dPres = fopen("data/dp","w");

    // ファイルが正しく開けたかを確認
    if (Ures == NULL) {
        printf("ファイルを開くことができません。\n");
        return; 
    }

    // ファイルに文字列を書き込む
    // U
    fprintf(Ures,"TimeStep=%d\n",timeStep);
    fprintf(dUres,"TimeStep=%d\n",timeStep);
    for(int i=0; i<yn+2; i++)
    {
        for(int j=0; j<xn+3; j++)
        {
            fprintf(Ures,"%.8e ",u[j+i*(xn+3)]);
            fprintf(dUres,"%.8e ",du[j+i*(xn+3)]);
        }
        fprintf(Ures, "\n");
        fprintf(dUres, "\n");
    }

    // V
    fprintf(Vres,"TimeStep=%d\n",timeStep);
    fprintf(dVres,"TimeStep=%d\n",timeStep);
    for(int i=0; i<yn+3; i++)
    {
        for(int j=0; j<xn+2; j++)
        {
            fprintf(Vres,"%.8e ",v[j+i*(xn+2)]);
            fprintf(dVres,"%.8e ",dv[j+i*(xn+2)]);
        }
        fprintf(Vres, "\n");
        fprintf(dVres, "\n");
    }
    // p
    fprintf(Pres,"TimeStep=%d\n",timeStep);
    fprintf(dPres,"TimeStep=%d\n",timeStep);
    for(int i=0; i<yn+2; i++)
    {
        for(int j=0; j<xn+2; j++)
        {
            // fprintf(Ures,"%.4lf ",(u[j+i*(xn+3)] + u[(j+1)+i*(xn+3)]) / 2.0);
            // fprintf(Vres,"%.4lf ",(v[j+i*(xn+2)] + v[j+(i+1)*(xn+2)]) / 2.0);
            fprintf(Pres,"%.8e ",p[j+i*(xn+2)]);
            fprintf(dPres,"%.8e ",dp[j+i*(xn+2)]);
        }
        // fprintf(Ures, "\n");
        // fprintf(Vres, "\n");
        fprintf(Pres, "\n");
        fprintf(dPres, "\n");
    }

    // ファイルを閉じる
    fclose(Ures);
    fclose(Vres);
    fclose(Pres);
    fclose(dUres);
    fclose(dVres);
    fclose(dPres);
}

// 残差などのログを出力する．
void writeLog(char* fileName, int timeStep, double dt, time_t elapsedTime, double* res, int resNum)
{
    // ログファイルを開く
    FILE* outputLog = fopen(fileName,"a");

    // ファイルが正しく開けたかを確認
    if (outputLog == NULL) {
        printf("Failed to open file\n");
        return; 
    }

    // 残差を書き込み
    // timeStep, time を書き込み
    fprintf(outputLog,"%d %.4lf",timeStep, timeStep*dt);

    // 残差を書き込み
    for (int i=0; i<resNum; i++)
    {
        fprintf(outputLog," %.8e",res[i]);
    }

    // 経過時間を書き込み
    fprintf(outputLog," %ld",elapsedTime);

    // 改行文字を書き込み
    fprintf(outputLog,"\n");

    // ファイルを閉じる
    fclose(outputLog);
}

// 残差を計算
// field        残差を計算する場(速度，圧力など)
// startIndexX,Y 残差計算に使用する領域の最初のインデックス．dUなら，X=2 (ゴーストセル，境界セルを除く) Y=1 (ゴーストセル除く)
// numX, numY   計算に使う領域の格子点数．基本はゴーストセルを除いた値を入れる．
// maxNumX      元の配列のx方向の要素数．dUならxn+3
// resStartIndex   残差をまとめて配列resに入れる．各フィールドがどの番号から始まるのか指定．Uは0,1, Vは2,3,... など．
void calcResidual(double* field, int startIndexX, int startIndexY, int numX, int numY, int maxNumX, double* res, int resStartIndex)
{
    double tmpErr = 0;
    for(int i=startIndexY; i<numY; i++)
    {
        for(int j=startIndexX; j<numX; j++)
        {   
            tmpErr += fabs(field[j+i*maxNumX]);
            // tmpErr += du[j+i*(numX+2)]*du[j+i*(numX+2)];
        }
    }
    
    // 残差を格納
    res[resStartIndex+0] = tmpErr;
    res[resStartIndex+1] = tmpErr/(numX*numY);
}



