#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265358979
#include <sys/stat.h> 

typedef struct _bmp 
{
    int Hpixels;                // 原圖像寬度（像素數）
    int Vpixels;                // 原圖像高度（像素數）
    unsigned char HeaderInfo[54];   // BMP Header (前 54 bytes)
    unsigned long int Hbytes;       // 一列(含padding)的 bytes 數
    unsigned char **data;           // 原始 BGR 資料(含padding)，大小為 [height][RowBytes]
    unsigned char **R;             // [height][width]
    unsigned char **G;             // [height][width]
    unsigned char **B;             // [height][width]
} bmp;

// 讀取 BMP 檔

int bmp_load_fn(const char *filename, bmp *p_bmp) 
{
    FILE* f = fopen(filename, "rb"); // 以 rb 打開 BMP
    if(!f){
        printf("\n\n%s NOT FOUND\n\n", filename);
        return 1;
    }
    fread(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f); // 讀取 BMP header 54 bytes
    
    int width, height;
    memcpy(&width,  &p_bmp->HeaderInfo[18], sizeof(int)); // offset=18
    memcpy(&height, &p_bmp->HeaderInfo[22], sizeof(int)); // offset=22
    p_bmp->Hpixels = width;
    p_bmp->Vpixels = height;

    // 每行實際的 byte 數 (要對齊到 4 bytes 邊界)
    int RowBytes = (width * 3 + 3) & (~3);
    p_bmp->Hbytes = RowBytes;

    // 配置記憶體給 data[]
    p_bmp->data = (unsigned char **)malloc(height * sizeof(unsigned char*));
    for(int i = 0; i < height; i++){
        p_bmp->data[i] = (unsigned char *)malloc(RowBytes * sizeof(unsigned char));
    }

    // 讀取每行像素資料
    for(int i = 0; i < height; i++){
        fread(p_bmp->data[i], sizeof(unsigned char), RowBytes, f);
    }
    fclose(f);

    // 配置 R, G, B 三個 2D array (大小就是原圖寬高)
    p_bmp->R = (unsigned char **)malloc(height * sizeof(unsigned char*));
    p_bmp->G = (unsigned char **)malloc(height * sizeof(unsigned char*));
    p_bmp->B = (unsigned char **)malloc(height * sizeof(unsigned char*));
    for(int i = 0; i < height; i++){
        p_bmp->R[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        p_bmp->G[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        p_bmp->B[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
    }

    // 分離出 R, G, B
    // BMP 預設像素順序為 B, G, R
    for(int i = 0; i < height; i++){
        for(int j = 0; j < width; j++){
            p_bmp->B[i][j] = p_bmp->data[i][3*j + 0];
            p_bmp->G[i][j] = p_bmp->data[i][3*j + 1];
            p_bmp->R[i][j] = p_bmp->data[i][3*j + 2];
        }
    }
    return 0; // success
}

// 釋放 bmp 的動態配置記憶體
int bmp_free(bmp *p_bmp)
{
    int height = p_bmp->Vpixels;
    // free data[]
    for(int i = 0; i < height; i++){
        free(p_bmp->data[i]);
        free(p_bmp->R[i]);
        free(p_bmp->G[i]);
        free(p_bmp->B[i]);
    }
    free(p_bmp->data);
    free(p_bmp->R);
    free(p_bmp->G);
    free(p_bmp->B);

    // 重置
    p_bmp->data = NULL;
    p_bmp->R = NULL;
    p_bmp->G = NULL;
    p_bmp->B = NULL;
    p_bmp->Hpixels = 0;
    p_bmp->Vpixels = 0;
    p_bmp->HeaderInfo[0] = '\0';
    p_bmp->Hbytes = 0;
    return 0;
}

// 配置 2D float, short
//=====================================================================
float **alloc2D_float(int h, int w)
{
    float **ptr = (float**)malloc(h * sizeof(float*));
    for(int i=0; i<h; i++){
        ptr[i] = (float*)malloc(w * sizeof(float));
    }
    return ptr;
}
void free2D_float(float **p, int h)
{
    for(int i=0; i<h; i++){
        free(p[i]);
    }
    free(p);
}
short **alloc2D_short(int h, int w)
{
    short **ptr = (short**)malloc(h * sizeof(short*));
    for(int i=0; i<h; i++){
        ptr[i] = (short*)malloc(w * sizeof(short));
    }
    return ptr;
}
void free2D_short(short **p, int h)
{
    for(int i=0; i<h; i++){
        free(p[i]);
    }
    free(p);
}

// 補邊: 若寬或高不是 8 的倍數, 補到 8 的倍數
// 這裡示範: Y, Cb, Cr 都會補 0
int getPaddedSize(int x) {    // 回傳最小的 >= x 的 8 倍數
    return ((x + 7) / 8) * 8;
}

// 將 (height, width) 大小的 channel 資料複製到 (padH, padW) 大小的新陣列
// 多餘的部分(右邊 & 下邊) 以 0 補
void copy_and_pad_float(float **dst, int padH, int padW,
                        float **src, int srcH, int srcW)
{
    for(int i=0; i<padH; i++){
        for(int j=0; j<padW; j++){
            if(i < srcH && j < srcW){
                dst[i][j] = src[i][j];
            }
            else {
                dst[i][j] = 0.0f; // padding
            }
        }
    }
}

// RGB -> YCbCr 
void rgb_to_ycbcr(bmp *p_bmp, float **Y, float **Cb, float **Cr)
{
    int h = p_bmp->Vpixels;
    int w = p_bmp->Hpixels;
    for(int i = 0; i < h; i++){
        for(int j = 0; j < w; j++){
            float R = (float)p_bmp->R[i][j];
            float G = (float)p_bmp->G[i][j];
            float B = (float)p_bmp->B[i][j];
            Y[i][j]  =  0.2126f * R + 0.7152f * G + 0.0722f * B;
            Cb[i][j] = -0.11457f * R - 0.38543f * G + 0.50000f * B + 128.0f;
            Cr[i][j] =  0.50000f * R - 0.45415f * G - 0.04585f * B + 128.0f;
        }
    }
}
// 計算 8x8 2D-DCT (簡易實作，非最佳化)
void dct_8x8(float inBlock[8][8], float outBlock[8][8])
{
    const float c0 = 1.0f / (float)sqrt(8.0);
    const float c1 = (float)sqrt(2.0/8.0);

    for(int u = 0; u < 8; u++){
        for(int v = 0; v < 8; v++){
            float cu = (u == 0) ? c0 : c1;
            float cv = (v == 0) ? c0 : c1;
            float sum = 0.0f;
            for(int m=0; m<8; m++){
                for(int n=0; n<8; n++){
                    float val = inBlock[m][n];
                    sum += val
                        * cos(((2*m + 1)*u*PI)/16.0f)
                        * cos(((2*n + 1)*v*PI)/16.0f);
                }
            }
            outBlock[u][v] = cu * cv * sum;
        }
    }
}
// 將整張 channel 的資料分成 8x8 block 做 DCT
// h, w 若是 8 的倍數 (這裡已經padding後)，即可安全覆蓋
void do_dct_for_channel(float **channel, float **Freq, int h, int w)
{
    int mb = h / 8;  // 垂直方向 block 數
    int nb = w / 8;  // 水平方向 block 數

    float inBlock[8][8], outBlock[8][8];

    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            // 收集 8x8
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 8; j++){
                    inBlock[i][j] = channel[ by*8 + i ][ bx*8 + j ];
                }
            }
            // 做 2D-DCT
            dct_8x8(inBlock, outBlock);
            // 放回 Freq
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 8; j++){
                    Freq[ by*8 + i ][ bx*8 + j ] = outBlock[i][j];
                }
            }
        }
    }
}
static int qtableY[8][8] = {
    {16, 11, 10, 16, 24,  40,  51,  61},
    {12, 12, 14, 19, 26,  58,  60,  55},
    {14, 13, 16, 24, 40,  57,  69,  56},
    {14, 17, 22, 29, 51,  87,  80,  62},
    {18, 22, 37, 56, 68, 109, 103,  77},
    {24, 35, 55, 64, 81, 104, 113,  92},
    {49, 64, 78, 87,103, 121, 120, 101},
    {72, 92, 95, 98,112, 100, 103,  99},
};

static int qtableCb[8][8] = {
    {17, 18, 24, 47, 99,  99,  99,  99},
    {18, 21, 26, 66, 99,  99,  99,  99},
    {24, 26, 56, 99, 99,  99,  99,  99},
    {47, 66, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
};
static int qtableCr[8][8] = {
    {17, 18, 24, 47, 99,  99,  99,  99},
    {18, 21, 26, 66, 99,  99,  99,  99},
    {24, 26, 56, 99, 99,  99,  99,  99},
    {47, 66, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
    {99, 99, 99, 99, 99,  99,  99,  99},
};

void quant_channel1(float **F, short **qF, float **eF,
                   int h, int w,
                   const int qtable[8][8])
{
    int mb = h / 8; 
    int nb = w / 8; 

    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    int row = by*8 + i;
                    int col = bx*8 + j;
                    float val = F[row][col];
                    int Q = qtable[i][j];
                    short qVal = (short)lrintf(val / (float)Q); // 四捨五入
                    qF[row][col] = qVal;

                    float recon = (float)qVal * (float)Q;
                    eF[row][col] = val - recon; // 誤差
                }
            }
        }
    }
}


// 量化 (F / Qt) -> 四捨五入存成 short
void quant_channel2(float **F, short **qF, int h, int w,
                   const int qtable[8][8])
{
    int mb = h / 8; 
    int nb = w / 8; 
    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    int row = by*8 + i;
                    int col = bx*8 + j;
                    float val = F[row][col];
                    int Q = qtable[i][j];
                    short qVal = (short)lrintf(val / (float)Q); 
                    qF[row][col] = qVal;
                }
            }
        }
    }
}

static const int ZigZagOrder[64] = {
     0,  1,  5,  6, 14, 15, 27, 28,
     2,  4,  7, 13, 16, 26, 29, 42,
     3,  8, 12, 17, 25, 30, 41, 43,
     9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};

void do_dpcm_for_channel(short **data, int padH, int padW)
{
    int mb = padH / 8;
    int nb = padW / 8;

    short prevDC = 0;
    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            short DC_cur = data[by*8][bx*8];     // block左上角位置 = DC
            short DC_diff = DC_cur - prevDC;
            data[by*8][bx*8] = DC_diff;         // 寫回此 block 的 DC
            prevDC = DC_cur;                    // 更新
        }
    }
}
// 將一個 8x8 block 的 64 筆係數做 ZigZag
void zigzag_block(const short block_in[64], short block_out[64])
{
    for(int i=0; i<64; i++){
        block_out[i] = block_in[ ZigZagOrder[i] ];
    }
}


// 將 8x8 量化表以 ASCII 輸出

void write_quant_table(const char *fname, const int qtable[8][8])
{
    FILE *fp = fopen(fname, "w");
    if(!fp){
        fprintf(stderr, "Cannot open %s\n", fname);
        return;
    }
    for(int i=0; i<8; i++){
        for(int j=0; j<8; j++){
            fprintf(fp, "%d ", qtable[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}



// RLE：將 64 筆 zigzag 後的資料編碼成 skip/value 形式
static void rle_encode_block(const short block[64], FILE *fp) {
    int skip = 0; 
    for (int i = 0; i < 64; i++) {
        short val = block[i];
        if (val == 0) {
            skip++;
        } else {
            if (skip > 0) {
                // 輸出 $skip N
                fprintf(fp, "$skip %d ", skip);
                skip = 0;
            }
            // 輸出 $value X
            fprintf(fp, "$value %d ", val);
        }
    }
    if (skip > 0) {   
        fprintf(fp, "$skip %d ", skip);
    }
    // block 結束換行
    fprintf(fp, "\n");
}

// 將一個 64 筆 zigzag 後的 block 寫成二進位 RLE
// 回傳寫到檔案的「byte 數量」
size_t rle_encode_block_binary(const short block[64], FILE *fp)
{
    size_t bytesWritten = 0;
    int skip = 0;
    for(int i=0; i<64; i++){
        short val = block[i];
        if(val == 0){
            skip++;
        } else {
            // 先把累積的 skip 寫出
            if(skip > 0){
                // op=0 => skip, 接下來 2 bytes = skip 數量
                unsigned char op = 0;
                fwrite(&op, 1, 1, fp);         // 1 byte
                fwrite(&skip, 2, 1, fp);       // 2 bytes (short)
                bytesWritten += 3;
                skip = 0;
            }
            // 再寫 value
            {
                // op=1 => value, 接下來 2 bytes = val
                unsigned char op = 1;
                fwrite(&op, 1, 1, fp);         // 1 byte
                fwrite(&val, 2, 1, fp);        // 2 bytes (short)
                bytesWritten += 3;
            }
        }
    }
    // 收尾處理可能剩餘的 skip
    if(skip > 0){
        unsigned char op = 0;  // skip
        fwrite(&op, 1, 1, fp);
        fwrite(&skip, 2, 1, fp);
        bytesWritten += 3;
    }
    return bytesWritten;
}

void do_dpcm_zigzag_rle(short **Y, short **Cb, short **Cr, 
                        int padH, int padW, FILE *fp_out)
{
    int mb = padH / 8; 
    int nb = padW / 8;

    // 1) 先對三個 channel 分別做 DPCM(DC)
    do_dpcm_for_channel(Y,  padH, padW);
    do_dpcm_for_channel(Cb, padH, padW);
    do_dpcm_for_channel(Cr, padH, padW);

    // 2) Block by Block 輸出
    short block_in[64], block_z[64];
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            //=== Y ===
            {
                int idx=0;
                for(int i=0;i<8;i++){
                    for(int j=0;j<8;j++){
                        block_in[idx++] = Y[by*8 + i][bx*8 + j];
                    }
                }
                zigzag_block(block_in, block_z);
                fprintf(fp_out, "(%d,%d, Y) ", by, bx);
                rle_encode_block(block_z, fp_out);
            }
            //=== Cb ===
            {
                int idx=0;
                for(int i=0;i<8;i++){
                    for(int j=0;j<8;j++){
                        block_in[idx++] = Cb[by*8 + i][bx*8 + j];
                    }
                }
                zigzag_block(block_in, block_z);
                fprintf(fp_out, "(%d,%d, Cb) ", by, bx);
                rle_encode_block(block_z, fp_out);
            }
            //=== Cr ===
            {
                int idx=0;
                for(int i=0;i<8;i++){
                    for(int j=0;j<8;j++){
                        block_in[idx++] = Cr[by*8 + i][bx*8 + j];
                    }
                }
                zigzag_block(block_in, block_z);
                fprintf(fp_out, "(%d,%d, Cr) ", by, bx);
                rle_encode_block(block_z, fp_out);
            }
        }
    }
}

// 對 Y/Cb/Cr 做 DC-DPCM -> ZigZag -> RLE(二進位)
// 回傳總共寫了幾個 bytes   不需在檔案中存放 block 索引，因為decoder 看得懂就好。
size_t do_dpcm_zigzag_rle_binary(short **Y, short **Cb, short **Cr,
                                 int padH, int padW, FILE *fp)
{
    int mb = padH / 8;
    int nb = padW / 8;
    // 先對三個 channel 做 DPCM (DC)
    // (掃描順序: 自左而右, 自上而下)
    size_t totalBytes = 0;

    // DPCM
    {
        short prevDC = 0;
        for(int by=0; by<mb; by++){
            for(int bx=0; bx<nb; bx++){
                short DC_cur = Y[by*8][bx*8];
                short DC_diff = DC_cur - prevDC;
                Y[by*8][bx*8] = DC_diff;
                prevDC = DC_cur;
            }
        }
    }
    {
        short prevDC = 0;
        for(int by=0; by<mb; by++){
            for(int bx=0; bx<nb; bx++){
                short DC_cur = Cb[by*8][bx*8];
                short DC_diff = DC_cur - prevDC;
                Cb[by*8][bx*8] = DC_diff;
                prevDC = DC_cur;
            }
        }
    }
    {
        short prevDC = 0;
        for(int by=0; by<mb; by++){
            for(int bx=0; bx<nb; bx++){
                short DC_cur = Cr[by*8][bx*8];
                short DC_diff = DC_cur - prevDC;
                Cr[by*8][bx*8] = DC_diff;
                prevDC = DC_cur;
            }
        }
    }

    // ZigZag + RLE
    short block_in[64], block_z[64];

    // === Y ===
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            int idx = 0;
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    block_in[idx++] = Y[by*8 + i][bx*8 + j];
                }
            }
            // ZigZag
            for(int i=0;i<64;i++){
                block_z[i] = block_in[ZigZagOrder[i]];
            }
            // RLE -> 寫檔
            size_t n = rle_encode_block_binary(block_z, fp);
            totalBytes += n;
        }
    }
    // === Cb ===
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            int idx = 0;
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    block_in[idx++] = Cb[by*8 + i][bx*8 + j];
                }
            }
            // ZigZag
            for(int i=0;i<64;i++){
                block_z[i] = block_in[ZigZagOrder[i]];
            }
            // RLE -> 寫檔
            size_t n = rle_encode_block_binary(block_z, fp);
            totalBytes += n;
        }
    }
    // === Cr ===
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            int idx = 0;
            for(int i=0;i<8;i++){
                for(int j=0;j<8;j++){
                    block_in[idx++] = Cr[by*8 + i][bx*8 + j];
                }
            }
            // ZigZag
            for(int i=0;i<64;i++){
                block_z[i] = block_in[ZigZagOrder[i]];
            }
            // RLE -> 寫檔
            size_t n = rle_encode_block_binary(block_z, fp);
            totalBytes += n;
        }
    }

    return totalBytes;
}



void save_quantized_short(const char *fname, short **qData, int h, int w)
{
    FILE *fp = fopen(fname, "wb");
    if(!fp){
        fprintf(stderr, "Cannot open %s\n", fname);
        return;
    }
    for(int i=0; i<h; i++){
        fwrite(qData[i], sizeof(short), w, fp);
    }
    fclose(fp);
}
void save_quantized_error(const char *fname, float **eData, int h, int w)
{
    FILE *fp = fopen(fname, "wb");
    if(!fp){
        fprintf(stderr, "Cannot open %s\n", fname);
        return;
    }
    for(int i=0; i<h; i++){
        fwrite(eData[i], sizeof(float), w, fp);
    }
    fclose(fp);
}
// 計算 SQNR(dB) = 10 log10( sum(F^2) / sum(eF^2) )
// 這裡「每一個頻率位置」都算一次 => 共 64 個結果
void compute_sqnr_8x8(float **F, float **eF, int h, int w, double sqnr[64])
{
    int mb = h / 8; 
    int nb = w / 8;

    double sumF2[64];
    double sumE2[64];
    for(int k=0; k<64; k++){
        sumF2[k] = 0.0;
        sumE2[k] = 0.0;
    }

    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            int baseY = by*8;
            int baseX = bx*8;
            for(int u=0; u<8; u++){
                for(int v=0; v<8; v++){
                    int freqIndex = u*8 + v;
                    float valF = F[baseY + u][baseX + v];
                    float vale = eF[baseY + u][baseX + v];
                    sumF2[freqIndex] += (double)valF * (double)valF;
                    sumE2[freqIndex] += (double)vale * (double)vale;
                }
            }
        }
    }
    for(int k=0; k<64; k++){
        if(sumE2[k] < 1e-12){ 
            sqnr[k] = 999.99; 
        }
        else {
            sqnr[k] = 10.0 * log10( sumF2[k] / sumE2[k] );
        }
    }
}
//將一個 8x8 block 的 64 個資料以 Run-Length Encoding 格式編碼
// 修正的 RLE 編碼函數，支持輸出 skip 和 value

void do_dpcm_for_all_blocks(short **data, int padH, int padW)
{
    int mb = padH / 8;
    int nb = padW / 8;

    short prevDC = 0; 
    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            // data的 block[0] 即為 row=(by*8), col=(bx*8)
            short DC_cur = data[by*8][bx*8];
            short DC_diff = DC_cur - prevDC;
            data[by*8][bx*8] = DC_diff;  // 用差分取代
            prevDC = DC_cur;            // 更新前一個 DC
        }
    }
}

void do_dpcm_zigzag_rle_for_channel(short **channel, int padH, int padW,
                                    const char *chName, FILE *fp)
{
    // 先針對 DC 做 DPCM
    do_dpcm_for_all_blocks(channel, padH, padW);

    int mb = padH / 8;
    int nb = padW / 8;
    short block_in[64], block_z[64];

    for (int by = 0; by < mb; by++) {
        for (int bx = 0; bx < nb; bx++) {
            // 收集當前 8x8 block
            int idx = 0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    block_in[idx++] = channel[by*8 + i][bx*8 + j];
                }
            }

            // ZigZag
            zigzag_block(block_in, block_z);

            // RLE 編碼
            // 在輸出前，標記一下 (by, bx, chName)
            fprintf(fp, "(%d,%d,%s) ", by, bx, chName);
            rle_encode_block(block_z, fp);
        }
    }
}
void do_dpcm_zigzag_rle_for_all(short **Y, short **Cb, short **Cr,
                                int padH, int padW, FILE *fp)
{
    do_dpcm_zigzag_rle_for_channel(Y,  padH, padW,  "Y",  fp);
    do_dpcm_zigzag_rle_for_channel(Cb, padH, padW,  "Cb", fp);
    do_dpcm_zigzag_rle_for_channel(Cr, padH, padW,  "Cr", fp);
}

unsigned int get_bmp_file_size(bmp *p_bmp) //取得 BMP 檔案大小 (Byte)
{
    // BMP header 第 2~5 bytes 是 file size
    unsigned int fileSize = 0;
    memcpy(&fileSize, &p_bmp->HeaderInfo[2], 4);
    return fileSize;
}

void print_compression_ratio(bmp *p_bmp,   //計算壓縮率
                            size_t Y_bytes, size_t Cb_bytes, size_t Cr_bytes,
                            size_t totalBytes_bin)
{
    int width  = p_bmp->Hpixels;
    int height = p_bmp->Vpixels;

    // 1) 取 BMP 原始大小
    unsigned int bmpFileSize = get_bmp_file_size(p_bmp); // 由 header 取得

    // 2) 計算三個 channel 原本大小 (不考慮 padding)
    double origY  = (double) (width * height);
    double origCb = (double) (width * height);
    double origCr = (double) (width * height);

    // 3) 壓縮率: origSize / compressedSize
    double ratioY  = origY  / (double)Y_bytes;
    double ratioCb = origCb / (double)Cb_bytes;
    double ratioCr = origCr / (double)Cr_bytes;

    // 4) 整體壓縮率: (BMP檔大小) / (rle_code.bin大小)
    double ratioAll = (double) bmpFileSize / (double) totalBytes_bin;

    // 螢幕列印
    printf("[Compression Ratio]\n");
    printf("Y-channel  : %.4f\n", ratioY);
    printf("Cb-channel : %.4f\n", ratioCb);
    printf("Cr-channel : %.4f\n", ratioCr);
    printf("Overall (BMP vs BIN): %.4f\n", ratioAll);
}


int main(int argc, char **argv)
{

    int mode = atoi(argv[1]);
    //   encoder 0 <inputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>
    if(mode == 0)
    {
        char *fn_input_bmp = argv[2];  // Kimberly.bmp
        char *fn_r         = argv[3];  // R.txt
        char *fn_g         = argv[4];  // G.txt
        char *fn_b         = argv[5];  // B.txt
        char *fn_dim       = argv[6];  // dim.txt

        bmp myBmp;
        if(bmp_load_fn(fn_input_bmp, &myBmp) != 0){
            return 1;
        }
        int width = myBmp.Hpixels;
        int height= myBmp.Vpixels;

        // 1) 將原圖寬高寫到 dim.txt
        FILE *fdim = fopen(fn_dim, "w");
        if(!fdim){
            fprintf(stderr, "Cannot open %s\n", fn_dim);
            return 1;
        }
        fprintf(fdim, "%d %d\n", width, height);
        fclose(fdim);

        // 2) 將 R channel 寫到 R.txt 
        FILE *fr = fopen(fn_r, "w");
        if(!fr){
            fprintf(stderr, "Cannot open %s\n", fn_r);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                fprintf(fr, "%d ", myBmp.R[i][j]);
            }
            fprintf(fr, "\n");
        }
        fclose(fr);

        // 3) 將 G channel 寫到 G.txt
        FILE *fg = fopen(fn_g, "w");
        if(!fg){
            fprintf(stderr, "Cannot open %s\n", fn_g);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                fprintf(fg, "%d ", myBmp.G[i][j]);
            }
            fprintf(fg, "\n");
        }
        fclose(fg);

        // 4) 將 B channel 寫到 B.txt
        FILE *fb = fopen(fn_b, "w");
        if(!fb){
            fprintf(stderr, "Cannot open %s\n", fn_b);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                fprintf(fb, "%d ", myBmp.B[i][j]);
            }
            fprintf(fb, "\n");
        }
        fclose(fb);

        bmp_free(&myBmp);
        printf("encoder (mode 0) done.\n");
    }
    // mode 1
    //   encoder 1 <inputBMP> <Qt_Y.txt> <Qt_Cb.txt> <Qt_Cr.txt> <dim.txt> 
    //             <qF_Y.raw> <qF_Cb.raw> <qF_Cr.raw> <eF_Y.raw> <eF_Cb.raw> <eF_Cr.raw>
    else if(mode == 1)
    {
        if(argc < 12){
            fprintf(stderr, "Arguments not enough for mode=1\n");
            return 1;
        }
        char *fn_input_bmp = argv[2];   // Kimberly.bmp
        char *fn_qtY       = argv[3];   // Qt_Y.txt
        char *fn_qtCb      = argv[4];   // Qt_Cb.txt
        char *fn_qtCr      = argv[5];   // Qt_Cr.txt
        char *fn_dim       = argv[6];   // dim.txt
        char *fn_qFY       = argv[7];   // qF_Y.raw
        char *fn_qFCb      = argv[8];   // qF_Cb.raw
        char *fn_qFCr      = argv[9];   // qF_Cr.raw
        char *fn_eFY       = argv[10];  // eF_Y.raw
        char *fn_eFCb      = argv[11];  // eF_Cb.raw
        char *fn_eFCr      = argv[12];  // eF_Cr.raw

        bmp myBmp;
        if(bmp_load_fn(fn_input_bmp, &myBmp) != 0){
            return 1;
        }
        int width  = myBmp.Hpixels;   // 原圖寬
        int height = myBmp.Vpixels;   // 原圖高

        // ------ 計算補邊後的尺寸 ------
        int padW = getPaddedSize(width);
        int padH = getPaddedSize(height);

        // 1) 將量化表存成 ASCII
        write_quant_table(fn_qtY,  qtableY);
        write_quant_table(fn_qtCb, qtableCb);
        write_quant_table(fn_qtCr, qtableCr);

        // 2) 將 (原圖寬高, pad後寬高) 寫到 dim.txt
        //    方便後續解碼或做對比時參考
        FILE *fdim = fopen(fn_dim, "w");
        if(!fdim){
            fprintf(stderr, "Cannot open %s\n", fn_dim);
            return 1;
        }
        fprintf(fdim, "%d %d\n", width, height);
        fclose(fdim);

        // 3) 配置 (原圖大小) 的 Y,Cb,Cr 暫存
        float **tmpY  = alloc2D_float(height, width);
        float **tmpCb = alloc2D_float(height, width);
        float **tmpCr = alloc2D_float(height, width);

        // 4) 先把 RGB->YCbCr 轉到 tmp*
        rgb_to_ycbcr(&myBmp, tmpY, tmpCb, tmpCr);

        // 5) 接著配置 (補邊後大小) 的 Y, Cb, Cr
        float **Y  = alloc2D_float(padH, padW);
        float **Cb = alloc2D_float(padH, padW);
        float **Cr = alloc2D_float(padH, padW);

        // 6) 把原圖的數據 copy 到新的 (pad 大小) 並補 0
        copy_and_pad_float(Y,  padH, padW, tmpY,  height, width);
        copy_and_pad_float(Cb, padH, padW, tmpCb, height, width);
        copy_and_pad_float(Cr, padH, padW, tmpCr, height, width);

        // 釋放暫存
        free2D_float(tmpY,  height);
        free2D_float(tmpCb, height);
        free2D_float(tmpCr, height);

        // 7) 對 (padH x padW) 做 2D-DCT
        float **FY  = alloc2D_float(padH, padW);
        float **FCb = alloc2D_float(padH, padW);
        float **FCr = alloc2D_float(padH, padW);

        do_dct_for_channel(Y,  FY,  padH, padW);
        do_dct_for_channel(Cb, FCb, padH, padW);
        do_dct_for_channel(Cr, FCr, padH, padW);

        // 8) 量化 + 誤差
        short **qFY  = alloc2D_short(padH, padW);
        short **qFCb = alloc2D_short(padH, padW);
        short **qFCr = alloc2D_short(padH, padW);

        float **eFY  = alloc2D_float(padH, padW);
        float **eFCb = alloc2D_float(padH, padW);
        float **eFCr = alloc2D_float(padH, padW);

        quant_channel1(FY,  qFY,  eFY,  padH, padW, qtableY);
        quant_channel1(FCb, qFCb, eFCb, padH, padW, qtableCb);
        quant_channel1(FCr, qFCr, eFCr, padH, padW, qtableCr);

        // 9) 輸出 qF_*.raw  (short)   eF_*.raw (float)
        //    注意: 檔案裡存的資料長寬 = (padH, padW)
        save_quantized_short(fn_qFY,  qFY,  padH, padW);
        save_quantized_short(fn_qFCb, qFCb, padH, padW);
        save_quantized_short(fn_qFCr, qFCr, padH, padW);

        save_quantized_error(fn_eFY,  eFY,  padH, padW);
        save_quantized_error(fn_eFCb, eFCb, padH, padW);
        save_quantized_error(fn_eFCr, eFCr, padH, padW);

        // 10) 計算 SQNR(dB): 每個頻率 0..63, 一共三個 channel
        double sqnrY[64], sqnrCb[64], sqnrCr[64];
        compute_sqnr_8x8(FY,  eFY,  padH, padW, sqnrY);
        compute_sqnr_8x8(FCb, eFCb, padH, padW, sqnrCb);
        compute_sqnr_8x8(FCr, eFCr, padH, padW, sqnrCr);

        // 11) 螢幕列印 3 x 64
        printf("SQNR(dB) for Y:\n");
        for(int k=0; k<64; k++){
            printf("%.2f ", sqnrY[k]);
        }
        printf("\n\n");

        printf("SQNR(dB) for Cb:\n");
        for(int k=0; k<64; k++){
            printf("%.2f ", sqnrCb[k]);
        }
        printf("\n\n");

        printf("SQNR(dB) for Cr:\n");
        for(int k=0; k<64; k++){
            printf("%.2f ", sqnrCr[k]);
        }
        printf("\n\n");

        // 收尾釋放
        free2D_float(Y,   padH);
        free2D_float(Cb,  padH);
        free2D_float(Cr,  padH);
        free2D_float(FY,  padH);
        free2D_float(FCb, padH);
        free2D_float(FCr, padH);

        free2D_short(qFY,  padH);
        free2D_short(qFCb, padH);
        free2D_short(qFCr, padH);

        free2D_float(eFY,  padH);
        free2D_float(eFCb, padH);
        free2D_float(eFCr, padH);

        bmp_free(&myBmp);

        printf("encoder (mode 1) done.\n");
    }
    else if (mode == 2) 
    {
        if (argc < 5) {
            fprintf(stderr, "Arguments not enough for mode=2\n");
            fprintf(stderr, "Usage: encoder 2 <bmpfile> ascii <rle_output>\n");
            return 1;
        }
        // argv[2]: inputBmp, argv[3]: "ascii", argv[4]: rle_output
        char *fn_input_bmp = argv[2];
        char *fmt_str      = argv[3];  // "ascii" (此處可不特別用，但根據題目要求先保留)
        char *fn_rle_out   = argv[4];

        // 讀取 BMP
        bmp myBmp;
        if (bmp_load_fn(fn_input_bmp, &myBmp) != 0) {
            return 1;
        }
        int width  = myBmp.Hpixels;
        int height = myBmp.Vpixels;

        // 打開輸出檔
        FILE *fp_out = fopen(fn_rle_out, (strcmp(fmt_str,"binary")==0)? "wb":"w");
    if (!fp_out) {
        fprintf(stderr, "Cannot open %s\n", fn_rle_out);
        return 1;
    }

    if(strcmp(fmt_str,"ascii")==0) {
        // ASCII 版：第一行印文字
        fprintf(fp_out, "%d %d\n", width, height);
    }
    else {
        // Binary 版：先寫 2 個 int 表示 width, height
        fwrite(&width,  sizeof(int), 1, fp_out);
        fwrite(&height, sizeof(int), 1, fp_out);
    }

        // 2) RGB->YCbCr
        float **tmpY  = alloc2D_float(height, width);
        float **tmpCb = alloc2D_float(height, width);
        float **tmpCr = alloc2D_float(height, width);
        rgb_to_ycbcr(&myBmp, tmpY, tmpCb, tmpCr);

        // 3) 補邊
        int padW = getPaddedSize(width);
        int padH = getPaddedSize(height);
        float **Y  = alloc2D_float(padH, padW);
        float **Cb = alloc2D_float(padH, padW);
        float **Cr = alloc2D_float(padH, padW);

        copy_and_pad_float(Y,  padH, padW, tmpY,  height, width);
        copy_and_pad_float(Cb, padH, padW, tmpCb, height, width);
        copy_and_pad_float(Cr, padH, padW, tmpCr, height, width);

        free2D_float(tmpY,  height);
        free2D_float(tmpCb, height);
        free2D_float(tmpCr, height);

        // 4) DCT
        float **FY  = alloc2D_float(padH, padW);
        float **FCb = alloc2D_float(padH, padW);
        float **FCr = alloc2D_float(padH, padW);
        do_dct_for_channel(Y,  FY,  padH, padW);
        do_dct_for_channel(Cb, FCb, padH, padW);
        do_dct_for_channel(Cr, FCr, padH, padW);

        // 5) 量化
        short **qFY  = alloc2D_short(padH, padW);
        short **qFCb = alloc2D_short(padH, padW);
        short **qFCr = alloc2D_short(padH, padW);
        quant_channel2(FY,  qFY,  padH, padW, qtableY);
        quant_channel2(FCb, qFCb, padH, padW, qtableCb);
        quant_channel2(FCr, qFCr, padH, padW, qtableCr);
        if(strcmp(fmt_str, "ascii")==0){
            do_dpcm_zigzag_rle(qFY, qFCb, qFCr, padH, padW, fp_out);
        }
        else {
        // -- (B) Binary 版 --
        //   先把 Y/Cb/Cr 通通做 DC-DPCM + ZigZag + RLE
        //   回傳總 bytes => 用來算壓縮率
        size_t totalBytesY  = 0;  
        size_t totalBytesCb = 0; 
        size_t totalBytesCr = 0; 

        // 我打算分開統計，所以拆成三段做：
        // (1) Y
        {
            // DPCM
            short prevDC = 0;
            int mb = padH/8, nb = padW/8;
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    short DC_cur = qFY[by*8][bx*8];
                    short DC_diff = DC_cur - prevDC;
                    qFY[by*8][bx*8] = DC_diff;
                    prevDC = DC_cur;
                }
            }
            // ZigZag+RLE
            short block_in[64], block_z[64];
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    int idx=0;
                    for(int i=0;i<8;i++){
                        for(int j=0;j<8;j++){
                            block_in[idx++] = qFY[by*8 + i][bx*8 + j];
                        }
                    }
                    for(int i=0;i<64;i++){
                        block_z[i] = block_in[ZigZagOrder[i]];
                    }
                    size_t n = rle_encode_block_binary(block_z, fp_out);
                    totalBytesY += n;
                }
            }
        }
        // (2) Cb
        {
            short prevDC = 0;
            int mb = padH/8, nb = padW/8;
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    short DC_cur = qFCb[by*8][bx*8];
                    short DC_diff = DC_cur - prevDC;
                    qFCb[by*8][bx*8] = DC_diff;
                    prevDC = DC_cur;
                }
            }
            short block_in[64], block_z[64];
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    int idx=0;
                    for(int i=0;i<8;i++){
                        for(int j=0;j<8;j++){
                            block_in[idx++] = qFCb[by*8 + i][bx*8 + j];
                        }
                    }
                    for(int i=0;i<64;i++){
                        block_z[i] = block_in[ZigZagOrder[i]];
                    }
                    size_t n = rle_encode_block_binary(block_z, fp_out);
                    totalBytesCb += n;
                }
            }
        }
        // (3) Cr
        {
            short prevDC = 0;
            int mb = padH/8, nb = padW/8;
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    short DC_cur = qFCr[by*8][bx*8];
                    short DC_diff = DC_cur - prevDC;
                    qFCr[by*8][bx*8] = DC_diff;
                    prevDC = DC_cur;
                }
            }
            short block_in[64], block_z[64];
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    int idx=0;
                    for(int i=0;i<8;i++){
                        for(int j=0;j<8;j++){
                            block_in[idx++] = qFCr[by*8 + i][bx*8 + j];
                        }
                    }
                    for(int i=0;i<64;i++){
                        block_z[i] = block_in[ZigZagOrder[i]];
                    }
                    size_t n = rle_encode_block_binary(block_z, fp_out);
                    totalBytesCr += n;
                }
            }
        }

        // 結算「總 RLE bytes」：檔案頭 (width+height=8 bytes) + 三個 channel
        size_t totalBytes_bin = 8 + totalBytesY + totalBytesCb + totalBytesCr;

        // 螢幕印出壓縮率
        //   - Y/Cb/Cr: (width*height) / (channel_bytes)
        //   - overall: BMPfileSize / rle_code.bin
        {
            double origY  = (double)(width * height);
            double origCb = (double)(width * height);
            double origCr = (double)(width * height);

            double ratioY  = origY  / (double)totalBytesY;
            double ratioCb = origCb / (double)totalBytesCb;
            double ratioCr = origCr / (double)totalBytesCr;

            // 取得 BMP 整體檔案大小
            unsigned int bmpFileSize = get_bmp_file_size(&myBmp);
            double ratioAll = (double)bmpFileSize / (double)totalBytes_bin;

            printf("[Compression Ratio]\n");
            printf("Y-channel  : %.4f\n", ratioY);
            printf("Cb-channel : %.4f\n", ratioCb);
            printf("Cr-channel : %.4f\n", ratioCr);
            printf("Overall (BMP vs BIN): %.4f\n", ratioAll);
        }
    }
        // 收尾
        fclose(fp_out);
        free2D_float(Y,  padH);   free2D_float(Cb,  padH);   free2D_float(Cr,  padH);
        free2D_float(FY, padH);   free2D_float(FCb, padH);   free2D_float(FCr, padH);
        free2D_short(qFY,  padH); free2D_short(qFCb, padH);  free2D_short(qFCr, padH);
        bmp_free(&myBmp);

        printf("encoder (mode 2) done.");
    }

    else {
        fprintf(stderr, "Unknown mode!\n");
    }

    return 0;
}

