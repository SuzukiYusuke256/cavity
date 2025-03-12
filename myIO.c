#include <stdio.h>
#include <stdarg.h> 
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "myIO.h"

int readConfig(const char* configName, int num, double* configArray)
{
    FILE* fp;
    char buffer[1024]; // configの1列目を読み込むためのバッファ
    
    // ファイルを開く
    fp = fopen(configName, "r");
    if (fp == NULL) {
        fprintf(stderr, "read config\nError: Failed to open %s\n", configName);
        return -1;
    }
    
    for (int ii = 0; ii < num; ii++) {
        fscanf(fp, "%s", buffer); // configの設定項目 (読み飛ばす)
        if (fscanf(fp, "%lf", &configArray[ii]) != 1) {
            fprintf(stderr, "read config\nError: Failed to read config data. Location: row %d\n", ii);
            fclose(fp);
            return -1;
        }
    }
    
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
            if (fscanf(fp, "%lf", &dataArray[y*numX + x]) != 1) {
                fprintf(stderr, "エラー: データの読み込みに失敗しました。位置: (%d, %d)\n", x, y);
                fclose(fp);
                return -1;
            }
        }
    }
    
    // ファイルを閉じる
    fclose(fp);
    return 0;
}

// データを書き込む
int writeData(char* fileName, double* field, int numX, int numY, char* header)
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
