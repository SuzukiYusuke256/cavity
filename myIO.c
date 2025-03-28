#include <stdio.h>
#include <stdarg.h> 
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include <dirent.h>
#include <ctype.h> // for isdigit

#include "myIO.h"
#include "myConst.h"
#include "config.h"

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
    config->writePrecision       = atoi(valueStrings[CFG_WRITE_PREC]);
    
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
int readData (double* field, int numX, int numY, const char* caseName, int timeStep, char* fieldName) 
{
    FILE* fp;
    char fileName[1024];
    char cTimeStep[128]; // store time as string
    char buffer[1024]; // ヘッダー行を読み込むためのバッファ
    
    // construct file name
    strcpy(fileName,caseName); // fileName = caseName
    strcat(fileName,"/"); // fileName = caseName/
    sprintf(cTimeStep,"%d",timeStep);
    strcat(fileName,cTimeStep); // fileName = caseName/timeStep
    strcat(fileName,"/"); // fileName = caseName/timeStep/
    strcat(fileName,fieldName); // fileName = caseName/timeStep/fieldName

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
            fscanf(fp, "%lf", &field[y*numX + x]);
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
int writeData(double* field, int numX, int numY, const char* caseName, int timeStep, char* fieldName, int nx, int ny, int writePrec)
{
    char fileName[1024];
    char cTimeStep[128]; // store time as string
    char buffer[1024]; // ヘッダー行を読み込むためのバッファ
    
    // construct file name
    strcpy(fileName,caseName);  // fileName = caseName
    strcat(fileName,"/");       // fileName = caseName/
    sprintf(cTimeStep,"%d",timeStep);
    strcat(fileName,cTimeStep); // fileName = caseName/timeStep
    strcat(fileName,"/");       // fileName = caseName/timeStep/
    strcat(fileName,fieldName); // fileName = caseName/timeStep/fieldName
    
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
            caseName,fieldName,timeStep,numX,numY,nx,ny);

    for(int jj=0; jj<numY; jj++)
    {
        for(int ii=0; ii<numX; ii++)
        {
            fprintf(fp,"%.*e ",writePrec,field[ii + jj*numX]); // writePrec: number of digits to write, field: data to write
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


/**
 * Retrieves a list of directories with numeric names within a specified case directory.
 * Converts directory names to double values and finds the maximum value.
 * 
 * @param caseDirPath       Path to the case directory
 * @param dirList           Array to store the numeric directory names as double values (must be pre-allocated)
 * @param maxDirs           Maximum number of directories to store in dirList
 * @param maxValue          Pointer to store the maximum numeric value found
 * @return Number of numeric directories found, or -1 on error
 */
int getNumericDirectories(const char* caseDirPath, double dirList[], int maxDirs, double* maxValue) 
{
    DIR* dir;
    struct dirent* entry;
    int count = 0;
    int tmpMaxValue = -1.0; // Initialize to a negative value to ensure any positive value will be greater

    // Open the directory
    dir = opendir(caseDirPath);
    if (dir == NULL) 
    {
        fprintf(stderr, "Error: Failed to open directory '%s'\n", caseDirPath);
        return -1;
    }

    // Iterate through directory entries
    while ((entry = readdir(dir)) != NULL) 
    {
        // Skip "." and ".."
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0) 
        {
            continue;
        }

        // Check if the name is numeric
        int isNumeric = 1;
        for (int i = 0; entry->d_name[i] != '\0'; i++) 
        {
            if (!isdigit(entry->d_name[i]) && entry->d_name[i] != '.') 
            {
                isNumeric = 0;
                break;
            }
        }

        // If numeric, convert to double and add to the list
        if (isNumeric && count < maxDirs) 
        {
            double value = atof(entry->d_name);
            dirList[count] = value;

            // Check if this is the maximum value so far
            if (value > tmpMaxValue) 
            {
                tmpMaxValue = value;
            }
            count++;
        }
    }

    // Set the maximum value and its corresponding directory name
    *maxValue = tmpMaxValue;

    closedir(dir);
    return count;
}