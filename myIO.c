#include <stdio.h>
#include <stdarg.h> 
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "myIO.h"
#include "myConst.h"

int readConfig(char* configName, Config* config)
{
    FILE* fp;
    // char buffer[1024]; // configの1列目を読み込むためのバッファ
    char line[1024];
    char *token;
    char *valueStrings[CFG_NUM];  // 文字列として保存する配列
    
    // ファイルを開く
    fp = fopen(configName, "r");
    if (fp == NULL) {
        fprintf(stderr, "read config\nError: Failed to open %s\n", configName);
        return -1;
    }

    int count = 0;
    // store setting as string
    while (fgets(line, 1024, fp) != NULL && count < CFG_NUM) {
        // 改行文字を削除
        line[strcspn(line, "\n")] = '\0';
        
        // 最初のトークン（パラメータ名）を取得
        token = strtok(line, " \t");
        if (token == NULL) continue; // 空行をスキップ
        
        // 2番目のトークン（値）を取得
        token = strtok(NULL, " \t");
        if (token == NULL) continue; // 2番目のトークンがない行をスキップ
        
        // 値を文字列配列にコピー
        valueStrings[count] = strdup(token);  // メモリ割り当てをしてコピー

        count++;
    }

    // copy to the struct
    strcpy(config->caseName,valueStrings[CFG_CASENAME]);
    config->stepNum              = atoi(valueStrings[CFG_STEP_NUM]);
    config->outputInterval       = atoi(valueStrings[CFG_OUTPUT_INT]);
    config->deltaT               = atof(valueStrings[CFG_DELTA_T]);
    config->nx                   = atoi(valueStrings[CFG_NX]);
    config->ny                   = atoi(valueStrings[CFG_NY]);
    config->re                   = atof(valueStrings[CFG_RE]);
    config->convergenceThreshold = atof(valueStrings[CFG_CONV_THRESH]);
    config->withInitialCondition = atoi(valueStrings[CFG_WITH_INIT]);
    
    // ファイルを閉じる
    fclose(fp);
    return 0;
}

/**
 * ファイルからデータを読み込む関数
 * 
 * @param fileName   読み込むファイルの名前
 * @param numX       x方向のデータ数
 * @param numY       y方向のデータ数
 * @param dataArray  読み込んだデータを格納する配列（サイズはnumX*numY）
 * @return           成功した場合は0、失敗した場合は-1を返す
 */
int readData(const char* fileName, int numX, int numY, double* dataArray) {
    FILE* fp;
    char buffer[1024]; // ヘッダー行を読み込むためのバッファ
    
    // ファイルを開く
    fp = fopen(fileName, "r");
    if (fp == NULL) {
        fprintf(stderr, "エラー: ファイル %s を開けませんでした。\n", fileName);
        return -1;
    }
    
    // ヘッダー行を読み飛ばす
    if (fgets(buffer, sizeof(buffer), fp) == NULL) {
        fprintf(stderr, "エラー: ヘッダー行を読み込めませんでした。\n");
        fclose(fp);
        return -1;
    }
    
    // データを読み込む
    for (int y = 0; y < numY; y++) {
        for (int x = 0; x < numX; x++) {
            fscanf(fp, "%lf", &dataArray[y*numX + x]);
            // printf("%ld\n",dataArray[y*numX + x]);

            // if (fscanf(fp, "%lf", &dataArray[y*numX + x]) != 1) {
            //     fprintf(stderr, "read Data\nError : データの読み込みに失敗しました。位置: (%d, %d)\n", x, y);
            //     fclose(fp);
            //     return -1;
            // }
        }
    }
    
    // ファイルを閉じる
    fclose(fp);
    return 0;
}

// データを書き込む
int writeDataHeader(char* fileName, double* field, int numX, int numY, char* header)
{
    // 出力
    FILE* fp = fopen(fileName,"w");

    // ファイルが正しく開けたかを確認
    if (fp == NULL) {
        printf("Failed to open file\n");
        return -1;
    }

    // ファイルに文字列を書き込む
    fprintf(fp, "%s\n",header);

    for(int jj=0; jj<numY; jj++)
    {
        for(int ii=0; ii<numX; ii++)
        {
            fprintf(fp,"%.8e ",field[ii + jj*numX]);
            // printf("%d %d %lf\n",ii,jj,field[ii + jj*numX]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}

// データを書き込む
// arrayNx, arrayNy : 
// nx, ny : number of cells in the computational domain.
int writeData(char* fileName, double* field, int arrayNx, int arrayNy, 
              char* caseName, char* fieldName, int nx, int ny, int timeStep)
{
    // 出力
    FILE* fp = fopen(fileName,"w");

    // ファイルが正しく開けたかを確認
    if (fp == NULL) {
        printf("Failed to open file\n");
        return -1;
    }

    // ファイルに文字列を書き込む
    char tmpStr[1024] = ""; 
    // sprintf(tmpStr,"%s sDUdy %dx%d (Nx:%d Ny:%d)",caseName,xn,1,xn,yn);
    fprintf(fp, "%s %s TimeStep:%d Size:%dx%d (Nx:%d Ny:%d)\n",
            caseName,fieldName,timeStep,arrayNx,arrayNy,nx,ny);

    for(int jj=0; jj<arrayNy; jj++)
    {
        for(int ii=0; ii<arrayNx; ii++)
        {
            fprintf(fp,"%.8e ",field[ii + jj*arrayNx]);
            // printf("%d %d %lf\n",ii,jj,field[ii + jj*numX]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    return 0;
}

/**
 * Creates a directory if it doesn't already exist
 * 
 * @param dirPath Path of the directory to create
 * @return 0 on success (directory exists or was created), -1 on failure
 */
int createDirectoryIfNotExists(const char* dirPath) {
    struct stat st;
    
    // Check if directory already exists
    if (stat(dirPath, &st) == 0) {
        if (S_ISDIR(st.st_mode)) {
            // Directory already exists
            return 0;
        } else {
            // Path exists but is not a directory
            fprintf(stderr, "Error: '%s' exists but is not a directory\n", dirPath);
            return -1;
        }
    }
    
    // Directory doesn't exist, create it
#ifdef _WIN32
    // Windows
    if (mkdir(dirPath) != 0) {
#else
    // Unix/Linux (permission 0755)
    if (mkdir(dirPath, 0755) != 0) {
#endif
        fprintf(stderr, "Error: Failed to create directory '%s'. Error code: %d\n", dirPath, errno);
        return -1;
    }
    
    printf("Created directory '%s'\n", dirPath);
    return 0;
}
