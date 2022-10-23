#include <mpi.h>
#include <stdio.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

#include "bmp.h"

using namespace std;

//定義平滑運算的次數
#define NSmooth 1000

/*********************************************************/
/*變數宣告：                                             */
/*  bmpHeader    ： BMP檔的標頭                          */
/*  bmpInfo      ： BMP檔的資訊                          */
/*  **BMPSaveData： 儲存要被寫入的像素資料               */
/*  **BMPData    ： 暫時儲存要被寫入的像素資料           */
/*********************************************************/
BMPHEADER bmpHeader;
BMPINFO bmpInfo;
RGBTRIPLE **BMPSaveData = NULL;
RGBTRIPLE **BMPData = NULL;

// Information for defining mpi type: mpi_rgbTriple
MPI_Datatype mpi_rgbTriple, mpi_rgbTripleRow;
constexpr const int BLOCK_LENGTH[3] = {1, 1, 1};
constexpr const MPI_Aint BLOCK_DISPLACEMENTS[3] = {0, 1, 2};
constexpr const MPI_Datatype BLOCK_TYPES[3] = {MPI_BYTE, MPI_BYTE, MPI_BYTE};

// indexer for bmpShape
constexpr int H = 0, W = 1;

/*********************************************************/
/*函數宣告：                                             */
/*  readBMP    ： 讀取圖檔，並把像素資料儲存在BMPSaveData*/
/*  saveBMP    ： 寫入圖檔，並把像素資料BMPSaveData寫入  */
/*  swap       ： 交換二個指標                           */
/*  **alloc_memory： 動態分配一個Y * X矩陣               */
/*********************************************************/
int readBMP(const char *fileName);  // read file
int saveBMP(const char *fileName);  // save file
void swap(RGBTRIPLE *a, RGBTRIPLE *b);
RGBTRIPLE **alloc_memory(int Y, int X);  // allocate memory
void transferBoundaries(int id, int comm_size, RGBTRIPLE **block, int blockHeight);

int main(int argc, char *argv[]) {
    /*********************************************************/
    /*變數宣告：                                             */
    /*  *infileName  ： 讀取檔名                             */
    /*  *outfileName ： 寫入檔名                             */
    /*  startTime   ： 記錄開始時間                         */
    /*  endwtime     ： 記錄結束時間                         */
    /*********************************************************/
    const char *infileName = "input.bmp";
    const char *outfileName = "output.bmp";
    double startTime = 0.0, endwtime = 0.;
    int id, comm_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    //記錄開始時間
    startTime = MPI_Wtime();

    //讀取檔案
    MPI_Datatype mpi_bmpShape;
    MPI_Type_contiguous(2, MPI_INT32_T, &mpi_bmpShape);
    MPI_Type_commit(&mpi_bmpShape);

    int bmpShape[2];
    
    if (id == 0) {
        if (readBMP(infileName))
            cout << "Read file successfully!!" << endl;
        else
            cout << "Read file fails!!" << endl;

        // broadcast image's shape to other processes.
        bmpShape[H] = bmpInfo.biHeight;
        bmpShape[W] = bmpInfo.biWidth;
        MPI_Bcast(bmpShape, 2, MPI_INT32_T, 0, MPI_COMM_WORLD);

        //動態分配記憶體給暫存空間
        BMPData = alloc_memory(bmpInfo.biHeight, bmpInfo.biWidth);
    }
    else {
        MPI_Bcast(bmpShape, 2, MPI_INT32_T, 0, MPI_COMM_WORLD);
    }


    // create type mpi_rgnTripleRow
    MPI_Type_create_struct(3, BLOCK_LENGTH, BLOCK_DISPLACEMENTS, BLOCK_TYPES, &mpi_rgbTriple);
    MPI_Type_contiguous(bmpShape[W], mpi_rgbTriple, &mpi_rgbTripleRow);
    MPI_Type_commit(&mpi_rgbTripleRow);

    int *sendcounts = new int[comm_size];
    int *displacements = new int[comm_size];
    int d = 0;
    for (int i = 0; i < comm_size; i++) {
        sendcounts[i] = (bmpShape[H] / comm_size) + (int)(i < bmpShape[H] % comm_size);
        displacements[i] = d;
        d += sendcounts[i];
    }
    
    if (id != 0) {
        BMPData = alloc_memory(sendcounts[id] + 2, bmpShape[W]);
        BMPSaveData = alloc_memory(sendcounts[id] + 2, bmpShape[W]);
    }
    else swap(BMPSaveData, BMPData);
    // preserve one row for boundary
    MPI_Scatterv(
        *BMPData, sendcounts, displacements, mpi_rgbTripleRow,
        *(&BMPSaveData[1]), sendcounts[id], mpi_rgbTripleRow, 0, MPI_COMM_WORLD
    );

    //進行多次的平滑運算
    for (int count = 0; count < NSmooth; count++) {

        transferBoundaries(id, comm_size, BMPSaveData, sendcounts[id] + 2);
        //把像素資料與暫存指標做交換
        swap(BMPSaveData, BMPData);
        //進行平滑運算
        for (int i = 1; i < sendcounts[id] + 1; i++)
            for (int j = 0; j < bmpShape[W]; j++) {
                /*********************************************************/
                /*設定上下左右像素的位置                                 */
                /*********************************************************/
                int Top = i > 0 ? i - 1 : sendcounts[id] + 2 - 1;
                int Down = i < sendcounts[id] + 2 - 1 ? i + 1 : 0;
                int Left = j > 0 ? j - 1 : bmpShape[W] - 1;
                int Right = j < bmpShape[W] - 1 ? j + 1 : 0;
                /*********************************************************/
                /*與上下左右像素做平均，並四捨五入                       */
                /*********************************************************/
                BMPSaveData[i][j].rgbBlue = (double)(BMPData[i][j].rgbBlue + BMPData[Top][j].rgbBlue + BMPData[Top][Left].rgbBlue + BMPData[Top][Right].rgbBlue + BMPData[Down][j].rgbBlue + BMPData[Down][Left].rgbBlue + BMPData[Down][Right].rgbBlue + BMPData[i][Left].rgbBlue + BMPData[i][Right].rgbBlue) / 9 + 0.5;
                BMPSaveData[i][j].rgbGreen = (double)(BMPData[i][j].rgbGreen + BMPData[Top][j].rgbGreen + BMPData[Top][Left].rgbGreen + BMPData[Top][Right].rgbGreen + BMPData[Down][j].rgbGreen + BMPData[Down][Left].rgbGreen + BMPData[Down][Right].rgbGreen + BMPData[i][Left].rgbGreen + BMPData[i][Right].rgbGreen) / 9 + 0.5;
                BMPSaveData[i][j].rgbRed = (double)(BMPData[i][j].rgbRed + BMPData[Top][j].rgbRed + BMPData[Top][Left].rgbRed + BMPData[Top][Right].rgbRed + BMPData[Down][j].rgbRed + BMPData[Down][Left].rgbRed + BMPData[Down][Right].rgbRed + BMPData[i][Left].rgbRed + BMPData[i][Right].rgbRed) / 9 + 0.5;
            }
    }

    MPI_Gatherv(
        *(&BMPSaveData[1]), sendcounts[id], mpi_rgbTripleRow,
        *BMPData, sendcounts, displacements, mpi_rgbTripleRow, 0, MPI_COMM_WORLD
    );
    swap(BMPData, BMPSaveData);
    //寫入檔案
    if (id == 0) {
        if (saveBMP(outfileName))
        // char outfileName_[20];
        // sprintf(outfileName_, "out_%2d.bmp", id);
        // if (saveBMP(outfileName_))
            cout << "Save file successfully!!" << endl;
        else
            cout << "Save file fails!!" << endl;
    }


    //得到結束時間，並印出執行時間
    if (id == 0) {
        endwtime = MPI_Wtime();
        cout << "The execution time = " << endwtime - startTime << endl ;
    }

    delete[] sendcounts;
    delete[] displacements;
    MPI_Type_free(&mpi_rgbTriple);
    MPI_Type_free(&mpi_rgbTripleRow);
    MPI_Type_free(&mpi_bmpShape);
    free(BMPData[0]);
    free(BMPSaveData[0]);
    free(BMPData);
    free(BMPSaveData);
    MPI_Finalize();

    return 0;
}

/*********************************************************/
/* 讀取圖檔                                              */
/*********************************************************/
int readBMP(const char *fileName) {
    //建立輸入檔案物件
    ifstream bmpFile(fileName, ios::in | ios::binary);

    //檔案無法開啟
    if (!bmpFile) {
        cout << "It can't open file!!" << endl;
        return 0;
    }

    //讀取BMP圖檔的標頭資料
    bmpFile.read((char *)&bmpHeader, sizeof(BMPHEADER));

    //判決是否為BMP圖檔
    if (bmpHeader.bfType != 0x4d42) {
        cout << "This file is not .BMP!!" << endl;
        return 0;
    }

    //讀取BMP的資訊
    bmpFile.read((char *)&bmpInfo, sizeof(BMPINFO));

    //判斷位元深度是否為24 bits
    if (bmpInfo.biBitCount != 24) {
        cout << "The file is not 24 bits!!" << endl;
        return 0;
    }

    //修正圖片的寬度為4的倍數
    while (bmpInfo.biWidth % 4 != 0)
        bmpInfo.biWidth++;

    //動態分配記憶體
    BMPSaveData = alloc_memory(bmpInfo.biHeight, bmpInfo.biWidth);

    //讀取像素資料
    // for(int i = 0; i < bmpInfo.biHeight; i++)
    //	bmpFile.read( (char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE));
    bmpFile.read((char *)BMPSaveData[0], bmpInfo.biWidth * sizeof(RGBTRIPLE) * bmpInfo.biHeight);

    //關閉檔案
    bmpFile.close();

    return 1;
}
/*********************************************************/
/* 儲存圖檔                                              */
/*********************************************************/
int saveBMP(const char *fileName) {
    //判決是否為BMP圖檔
    if (bmpHeader.bfType != 0x4d42) {
        cout << "This file is not .BMP!!" << endl;
        return 0;
    }

    //建立輸出檔案物件
    ofstream newFile(fileName, ios::out | ios::binary);

    //檔案無法建立
    if (!newFile) {
        cout << "The File can't create!!" << endl;
        return 0;
    }

    //寫入BMP圖檔的標頭資料
    newFile.write((char *)&bmpHeader, sizeof(BMPHEADER));

    //寫入BMP的資訊
    newFile.write((char *)&bmpInfo, sizeof(BMPINFO));

    //寫入像素資料
    // for( int i = 0; i < bmpInfo.biHeight; i++ )
    //        newFile.write( ( char* )BMPSaveData[i], bmpInfo.biWidth*sizeof(RGBTRIPLE) );
    newFile.write((char *)BMPSaveData[0], bmpInfo.biWidth * sizeof(RGBTRIPLE) * bmpInfo.biHeight);

    //寫入檔案
    newFile.close();

    return 1;
}

/*********************************************************/
/* 分配記憶體：回傳為Y*X的矩陣                           */
/*********************************************************/
RGBTRIPLE **alloc_memory(int Y, int X) {
    //建立長度為Y的指標陣列
    RGBTRIPLE **temp = new RGBTRIPLE *[Y];
    RGBTRIPLE *temp2 = new RGBTRIPLE[Y * X];
    memset(temp, 0, sizeof(RGBTRIPLE) * Y);
    memset(temp2, 0, sizeof(RGBTRIPLE) * Y * X);

    //對每個指標陣列裡的指標宣告一個長度為X的陣列
    for (int i = 0; i < Y; i++) {
        temp[i] = &temp2[i * X];
    }

    return temp;
}
/*********************************************************/
/* 交換二個指標                                          */
/*********************************************************/
void swap(RGBTRIPLE *a, RGBTRIPLE *b) {
    RGBTRIPLE *temp;
    temp = a;
    a = b;
    b = temp;
}

void transferBoundaries(int id, int comm_size, RGBTRIPLE **block, int blockHeight) {

    int next_id = (id + 1) % comm_size;
    int prev_id = (id + comm_size - 1) % comm_size;

    MPI_Send(*(&block[1]), 1, mpi_rgbTripleRow, prev_id, 0, MPI_COMM_WORLD);
    MPI_Send(*(&block[blockHeight - 2]), 1, mpi_rgbTripleRow, next_id, 0, MPI_COMM_WORLD);
    MPI_Recv(*(&block[0]), 1, mpi_rgbTripleRow, prev_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(*(&block[blockHeight - 1]), 1, mpi_rgbTripleRow, next_id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

