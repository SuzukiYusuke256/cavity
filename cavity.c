// cavity flow solver using MAC method
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// 境界条件に合わせてセルの値を更新
void correctBoundaryConditions(double* u, double* v, double* p, int xn, int yn)
{
    // 上下面
    double uWall = 1.0;
    for(int i=1; i<xn+2; i++)
    {
        // 壁面速度による境界条件
        u[i + (yn+1) * (xn+3)] = 2.0 * uWall - u[i + yn * (xn+3)];
        u[i + 0 * (xn+3)] = 2.0 * 0 - u[i + 1 * (xn+3)]; // すべりなし
    }
    for(int i=1; i<xn+1; i++)
    {
        // 壁面速度による境界条件
        v[i + 1 * (xn+2)] = 0;
        v[i + (yn+1) * (xn+2)] = 0;

        // 連続の式による境界条件
        v[i + 0 * (xn+2)] = v[i + 2 * (xn+2)]; 
        v[i + (yn+2) * (xn+2)] = v[i + yn * (xn+2)];
    }

    // 左右面
    for(int j=1; j<yn+2; j++)
    {
        // 壁面速度による境界条件
        v[0 + j * (xn+2)] = 2.0 * 0 - v[1 + j * (xn+2)]; // すべりなし
        v[(xn+1) + j * (xn+2)] = 2.0 * 0 - v[xn + j * (xn+2)]; // すべりなし
    }
    for(int j=1; j<xn+1; j++)
    {
        // 壁面速度による境界条件
        u[1 + j * (xn+3)] = 0;
        u[(xn+1) + j * (xn+3)] = 0; // すべりなし
    
        // 連続の式による境界条件
        u[0 + j * (xn+3)] = u[2 + j * (xn+3)];
        u[(xn+2) + j * (xn+3)] = u[xn + j * (xn+3)];
    }

    // 圧力境界条件
    // 壁面での勾配0
    // 上下面
    for(int i=0; i<xn+2; i++)
    {
        p[i + 0*(xn+2)] = p[i + 1*(xn+2)];
        p[i + (yn+1)*(xn+2)] = p[i + yn*(xn+2)];
    }
    //左右面
    for(int j=0; j<yn+2; j++)
    {
        p[0 + j * (xn+2)] = p[1 + j * (xn+2)];
        p[(xn+1) + j * (xn+2)] = p[xn + j * (xn+2)];
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
        err = 0.0;

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
void correct(double* u, double* v, double* p, double* du, double* dv, double* dp, double dx, double dy, double dt, int xn, int yn)
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
    for(int i=1; i<yn+1; i++)
    {
        for(int j=1; j<xn+1; j++)
        {
            p[j+i*(xn+2)] += dp[j+i*(xn+2)];
        }
    }
}

// データを書き込む
void write(double* u, double* v, double* p, double* du, double* dv, double* dp, int xn, int yn, char* dataDirName)
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
    for(int i=0; i<yn+2; i++)
    {
        for(int j=0; j<xn+3; j++)
        {
            fprintf(Ures,"%.4lf ",u[j+i*(xn+3)]);
            // fprintf(dUres,"%.4lf ",du[j+i*(xn+3)]);
        }
        fprintf(Ures, "\n");
        // fprintf(dUres, "\n");
    }

    fclose(Ures);
    return;

    // V
    for(int i=0; i<yn+3; i++)
    {
        for(int j=0; j<xn+2; j++)
        {
            fprintf(Vres,"%.4lf ",v[j+i*(xn+2)]);
            fprintf(dVres,"%.4lf ",dv[j+i*(xn+2)]);
        }
        fprintf(Vres, "\n");
        fprintf(dVres, "\n");
    }
    // p
    for(int i=0; i<yn+2; i++)
    {
        for(int j=0; j<xn+2; j++)
        {
            // fprintf(Ures,"%.4lf ",(u[j+i*(xn+3)] + u[(j+1)+i*(xn+3)]) / 2.0);
            // fprintf(Vres,"%.4lf ",(v[j+i*(xn+2)] + v[j+(i+1)*(xn+2)]) / 2.0);
            fprintf(Pres,"%.4lf ",p[j+i*(xn+2)]);
            fprintf(dPres,"%.4lf ",dp[j+i*(xn+2)]);
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
void writeLog(char* fileName, int timeStep, double dt, double* res, int resNum)
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
        fprintf(outputLog," %.4e",res[i]);
    }
    // 改行文字を書き込み
    fprintf(outputLog,"\n");

    // ファイルを閉じる
    fclose(outputLog);
}

// 残差を計算
// field        残差を計算する場(速度，圧力など)
// numX, numY   計算に使う領域の格子点数．基本はゴーストセルを除いた値を入れる．
// startIndex   残差をまとめて配列resに入れる．各フィールドがどの番号から始まるのか指定．Uは0,1, Vは2,3,... など．
void calcResidual(double* field, int numX, int numY, double* res, int startIndex)
{
    double tmpErr = 0;
    for(int i=1; i<numY+1; i++)
    {
        for(int j=1; j<numX+1; j++)
        {   
            tmpErr += fabs(field[j+i*(numX+2)]);
            // tmpErr += du[j+i*(numX+2)]*du[j+i*(numX+2)];
        }
    }
    
    // 残差を格納
    res[startIndex+0] = tmpErr;
    res[startIndex+1] = tmpErr/(numX*numY);
}

int main()
{
    // 計算設定
    int stepNum = 1;
    int outputInterval = 10;

    double dt = 0.001;

    int xn = 128;
    int yn = 128;

    // 流体の物性値
    double Re = 100;

    // output file name
    char* logFileName = "data/log";

    double dx = 1.0 / (double)xn;
    double dy = 1.0 / (double)yn;
    
    double* u = (double*)calloc((xn+3)*(yn+2),sizeof(double));
    double* v = (double*)calloc((xn+2)*(yn+3),sizeof(double));
    double* p = (double*)calloc((xn+2)*(yn+2),sizeof(double));
    double* du = (double*)calloc((xn+3)*(yn+2),sizeof(double));
    double* dv = (double*)calloc((xn+2)*(yn+3),sizeof(double));
    double* dp = (double*)calloc((xn+2)*(yn+2),sizeof(double));

    double residuals[6] = {0}; // ログに残差を出力するための配列．0で初期化

    // 境界条件
    // y = yn+2 u 速度一定の壁面
    // y = 1, x=1,xn+2 速度0
    // x = 0,xn+3 y = 0,y+3 連続の式を満たすように設定

    // initial conditions
    // U
    for(int i=1; i<yn+1; i++)
    {
        for(int j=1; j<xn+2; j++)
        {
            u[j+i*(xn+3)] = 0;
        }
    }
    // V
    for(int i=1; i<yn+2; i++)
    {
        for(int j=1; j<xn+1; j++)
        {
            v[j+i*(xn+2)] = 0;
        }
    }

    // boundary conditions
    double time = 0;
    char timeDirName[100] = {0,}; // 使っていない

    // ログファイルのヘッダを書き込み
    FILE* outputLog = fopen(logFileName,"w");

    // ファイルが正しく開けたかを確認
    if (outputLog == NULL) {
        printf("Failed to open file\n");
        return 0;
    }
    else
    {
        fprintf(outputLog,"TimeStep Time URes URes(cell) VRes VRes(cell) pRes pRes(cell)\n"); // ヘッダを書き込み
        fclose(outputLog); // ファイルを閉じる
    }

    // SMAC法での計算
    for(int k=0; k<=stepNum; k++)
    {
        // 境界条件の設定
        correctBoundaryConditions(u,v,p,xn,yn);

        // 仮速度の計算
        // u,vに仮速度u*,v*が入る
        predictVelocity(u,v,p,du,dv,dp,Re,dx,dy,dt,xn,yn);

        // 圧力ポアソン方程式を解く
        // 圧力修正子を計算
        poissonPressure(u,v,p,du,dv,dp,Re,dx,dy,dt,xn,yn);

        // printf("called\n");

        // 速度，圧力の修正
        correct(u,v,p,du,dv,dp,dx,dy,dt,xn,yn);

        // printf("called\n");

        // write(u,v,p,du,dv,dp,xn,yn,timeDirName);

        // 残差を出力
        if (k % outputInterval == 0)
        {
            // 残差の計算
            calcResidual(du,xn+1,yn,residuals,0);
            calcResidual(dv,xn,yn+1,residuals,2);
            calcResidual(dp,xn,yn,residuals,4);

            // ログに残差を出力
            writeLog(logFileName,k,dt,residuals,6);
        }

        time = time + dt;

    }
    
    strcpy(timeDirName,"data");
    // printf("%s\n",timeDirName);

    // //出力
    // write(u,v,p,du,dv,dp,xn,yn,);

    // 関数の引数には char* をいれないとだめ？  
    write(u,v,p,du,dv,dp,xn,yn,timeDirName);

    // メモリの解放
    free(u);
    free(v);
    free(p);
    free(du);
    free(dv);
    free(dp);

    // printf("%lf\n",u[3]);

}

