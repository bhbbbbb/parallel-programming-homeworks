#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
// #include <ctime>
#include <chrono>
#include <pthread.h>
#include <semaphore.h>

#include "bmp.h"


//定義平滑運算的次數
#define NSmooth 1000
// #define NSmooth 10

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
int num_thread;
int counter;
sem_t barrier_sems[2];
pthread_mutex_t counter_mutex;


/*********************************************************/
/*函數宣告：                                             */
/*  readBMP    ： 讀取圖檔，並把像素資料儲存在BMPSaveData*/
/*  saveBMP    ： 寫入圖檔，並把像素資料BMPSaveData寫入  */
/*  swap       ： 交換二個指標                           */
/*  **alloc_memory： 動態分配一個Y * X矩陣               */
/*********************************************************/
int readBMP(const char *fileName);  // read file
int saveBMP(const char *fileName);  // save file
RGBTRIPLE **alloc_memory(int Y, int X);  // allocate memory
void* multi_smooth(void* arg);
inline void _smooth(int rank);

int main(int argc, char *argv[]) {
    /*********************************************************/
    /*變數宣告：                                             */
    /*  *infileName  ： 讀取檔名                             */
    /*  *outfileName ： 寫入檔名                             */
    /*  startwtime   ： 記錄開始時間                         */
    /*  endwtime     ： 記錄結束時間                         */
    /*********************************************************/
    const char *infileName = "input.bmp";
    const char *outfileName = "output.bmp";

    if (readBMP(infileName))
        std::cout << "Read file successfully!!" << std::endl;
    else
        std::cout << "Read file fails!!" << std::endl;

    auto start_time = std::chrono::system_clock::now();

    // *********** shared variables initialization ***************
    num_thread = argc < 2 ? 1 : std::atoi(argv[1]);
    BMPData = alloc_memory(bmpInfo.biHeight, bmpInfo.biWidth); // empty here
    std::swap(BMPData, BMPSaveData);
    std::vector<pthread_t> threads(num_thread);
    pthread_mutex_init(&counter_mutex, NULL);

    for (auto& barrier_sem : barrier_sems)
        sem_init(&barrier_sem, 1, 0);

    counter = 0;
    // *********** shared variables initialization ***************

    for (long t = 0; t < num_thread; t++)
        pthread_create(&threads[t], NULL, multi_smooth, (void*)t);
    
    for (auto& t : threads)
        pthread_join(t, NULL);
    
    //得到結束時間，並印出執行時間
    std::chrono::duration<double> duration = std::chrono::system_clock::now() - start_time;
    std::cout << "The execution time = " << duration.count() << std::endl ;
    std::printf("n = %d, time = %f s", num_thread, duration.count());
    
    //寫入檔案
    if (saveBMP(outfileName))
        std::cout << "Save file successfully!!" << std::endl;
    else
        std::cout << "Save file fails!!" << std::endl;


    free(BMPData[0]);
    free(BMPSaveData[0]);
    free(BMPData);
    free(BMPSaveData);
    pthread_mutex_destroy(&counter_mutex);
    for (auto& barrier_sem : barrier_sems)
        sem_destroy(&barrier_sem);
    return 0;
}

/*********************************************************/
/* 讀取圖檔                                              */
/*********************************************************/
int readBMP(const char *fileName) {
    //建立輸入檔案物件
    std::ifstream bmpFile(fileName, std::ios::in | std::ios::binary);

    //檔案無法開啟
    if (!bmpFile) {
        std::cout << "It can't open file!!" << std::endl;
        return 0;
    }

    //讀取BMP圖檔的標頭資料
    bmpFile.read((char *)&bmpHeader, sizeof(BMPHEADER));

    //判決是否為BMP圖檔
    if (bmpHeader.bfType != 0x4d42) {
        std::cout << "This file is not .BMP!!" << std::endl;
        return 0;
    }

    //讀取BMP的資訊
    bmpFile.read((char *)&bmpInfo, sizeof(BMPINFO));

    //判斷位元深度是否為24 bits
    if (bmpInfo.biBitCount != 24) {
        std::cout << "The file is not 24 bits!!" << std::endl;
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
        std::cout << "This file is not .BMP!!" << std::endl;
        return 0;
    }

    //建立輸出檔案物件
    std::ofstream newFile(fileName, std::ios::out | std::ios::binary);

    //檔案無法建立
    if (!newFile) {
        std::cout << "The File can't create!!" << std::endl;
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

inline void _smooth(long rank) {

    for (int i = rank; i < bmpInfo.biHeight; i += num_thread)
        for (int j = 0; j < bmpInfo.biWidth; j++) {
            /*********************************************************/
            /*設定上下左右像素的位置                                 */
            /*********************************************************/
            int Top = i > 0 ? i - 1 : bmpInfo.biHeight - 1;
            int Down = i < bmpInfo.biHeight - 1 ? i + 1 : 0;
            int Left = j > 0 ? j - 1 : bmpInfo.biWidth - 1;
            int Right = j < bmpInfo.biWidth - 1 ? j + 1 : 0;
            /*********************************************************/
            /*與上下左右像素做平均，並四捨五入                       */
            /*********************************************************/
            BMPSaveData[i][j].rgbBlue = (double)(BMPData[i][j].rgbBlue + BMPData[Top][j].rgbBlue + BMPData[Top][Left].rgbBlue + BMPData[Top][Right].rgbBlue + BMPData[Down][j].rgbBlue + BMPData[Down][Left].rgbBlue + BMPData[Down][Right].rgbBlue + BMPData[i][Left].rgbBlue + BMPData[i][Right].rgbBlue) / 9 + 0.5;
            BMPSaveData[i][j].rgbGreen = (double)(BMPData[i][j].rgbGreen + BMPData[Top][j].rgbGreen + BMPData[Top][Left].rgbGreen + BMPData[Top][Right].rgbGreen + BMPData[Down][j].rgbGreen + BMPData[Down][Left].rgbGreen + BMPData[Down][Right].rgbGreen + BMPData[i][Left].rgbGreen + BMPData[i][Right].rgbGreen) / 9 + 0.5;
            BMPSaveData[i][j].rgbRed = (double)(BMPData[i][j].rgbRed + BMPData[Top][j].rgbRed + BMPData[Top][Left].rgbRed + BMPData[Top][Right].rgbRed + BMPData[Down][j].rgbRed + BMPData[Down][Left].rgbRed + BMPData[Down][Right].rgbRed + BMPData[i][Left].rgbRed + BMPData[i][Right].rgbRed) / 9 + 0.5;
        }
}

void* multi_smooth(void* arg) {

    long rank = reinterpret_cast<long>(arg);

    //進行多次的平滑運算
    for (int count = 0; count < NSmooth; count++) {
        //把像素資料與暫存指標做交換

        _smooth(rank);

        pthread_mutex_lock(&counter_mutex);

        if (counter < num_thread - 1) {

            counter++;
            pthread_mutex_unlock(&counter_mutex);
            sem_wait(&barrier_sems[count & 1]);
        }

        else {
            // counter == num_thread - 1
            counter = 0;
            pthread_mutex_unlock(&counter_mutex);

            if (count < NSmooth - 1) std::swap(BMPSaveData, BMPData);
            for (int i = 0; i < num_thread - 1; i++)
                sem_post(&barrier_sems[count & 1]);
        }
    }

    return NULL;
}