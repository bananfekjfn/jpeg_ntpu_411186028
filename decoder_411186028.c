#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265358979
#include <math.h>

// 結構與通用函式
typedef struct _bmp
{
    int Hpixels; // 圖片寬度
    int Vpixels; 
    unsigned char HeaderInfo[54];
    unsigned long int Hbytes;
    unsigned char **data; // 二維陣列儲存像素資料 (含 padding)
    unsigned char **R;    // 二維陣列儲存 R 
    unsigned char **G; 
    unsigned char **B;
} bmp;

// 配置、釋放 2D unsigned char
unsigned char **alloc2D_uc(int h, int w)
{
    unsigned char **p = (unsigned char**)malloc(h * sizeof(unsigned char*));
    for(int i=0; i<h; i++){
        p[i] = (unsigned char*)malloc(w * sizeof(unsigned char));
    }
    return p;
}

void free2D_uc(unsigned char **p, int h)
{
    for(int i=0; i<h; i++){
        free(p[i]);
    }
    free(p);
}

// 配置、釋放 2D float
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

// 配置、釋放 2D short
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
// BMP Header 初始化、儲存與釋放
void init_bmp_header(bmp *p_bmp) 
{
    // 先全部清 0
    memset(p_bmp->HeaderInfo, 0, 54);

    // "BM"
    p_bmp->HeaderInfo[0] = 'B'; 
    p_bmp->HeaderInfo[1] = 'M';
    
    // 檔案總大小 = Header(54) + 像素資料
    int fileSize = (p_bmp->Hbytes * p_bmp->Vpixels) + 54;
    memcpy(&p_bmp->HeaderInfo[2], &fileSize, 4);

    // 像素資料起始位移
    int bfOffBits = 54; 
    memcpy(&p_bmp->HeaderInfo[10], &bfOffBits, 4);

    // InfoHeader 大小
    int biSize = 40; 
    memcpy(&p_bmp->HeaderInfo[14], &biSize, 4);

    // 寬高
    memcpy(&p_bmp->HeaderInfo[18], &p_bmp->Hpixels, 4);
    memcpy(&p_bmp->HeaderInfo[22], &p_bmp->Vpixels, 4);

    // planes = 1
    short biPlanes = 1;
    memcpy(&p_bmp->HeaderInfo[26], &biPlanes, 2);

    // bitCount = 24
    short biBitCount = 24;
    memcpy(&p_bmp->HeaderInfo[28], &biBitCount, 2);

    // 影像資料大小
    int biSizeImage = p_bmp->Hbytes * p_bmp->Vpixels;
    memcpy(&p_bmp->HeaderInfo[34], &biSizeImage, 4);

    // X/Y 方向 DPI (2835 ~ 72 DPI)
    int biXPelsPerMeter = 2835; 
    memcpy(&p_bmp->HeaderInfo[38], &biXPelsPerMeter, 4);

    int biYPelsPerMeter = 2835;
    memcpy(&p_bmp->HeaderInfo[42], &biYPelsPerMeter, 4);
}

int bmp_load_fn(char *filename, bmp *p_bmp)
{
    FILE *f = fopen(filename, "rb");
    if(!f){
        fprintf(stderr, "Cannot open BMP file: %s\n", filename);
        return 0;
    }
    // 先讀前 54 bytes
    fread(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f);

    int width, height;
    memcpy(&width,  &p_bmp->HeaderInfo[18], sizeof(int)); // offset=18
    memcpy(&height, &p_bmp->HeaderInfo[22], sizeof(int)); // offset=22
    p_bmp->Hpixels = width;
    p_bmp->Vpixels = height;

    // 每行實際的 byte 數 (要對齊到 4 bytes 邊界)
    int RowBytes = (width * 3 + 3) & (~3);
    p_bmp->Hbytes = RowBytes;

    // 分配記憶體
    p_bmp->data = alloc2D_uc(p_bmp->Vpixels, p_bmp->Hbytes);
    p_bmp->R    = alloc2D_uc(p_bmp->Vpixels, p_bmp->Hpixels);
    p_bmp->G    = alloc2D_uc(p_bmp->Vpixels, p_bmp->Hpixels);
    p_bmp->B    = alloc2D_uc(p_bmp->Vpixels, p_bmp->Hpixels);

    // 讀像素資料
    for(int i=0; i < p_bmp->Vpixels; i++){
        fread(p_bmp->data[i], sizeof(unsigned char), p_bmp->Hbytes, f);
    }
    fclose(f);

    // 分離 BGR -> R/G/B
    for(int i=0; i < p_bmp->Vpixels; i++){
        for(int j=0; j < p_bmp->Hpixels; j++){
            p_bmp->B[i][j] = p_bmp->data[i][3*j + 0];
            p_bmp->G[i][j] = p_bmp->data[i][3*j + 1];
            p_bmp->R[i][j] = p_bmp->data[i][3*j + 2];
        }
    }
    return 1;
}


int bmp_save_fn(char *filename, bmp *p_bmp) 
{
    FILE* f = fopen(filename, "wb");
    if(f == NULL){
        printf("\n\nFILE CREATION ERROR: %s\n\n", filename);
        exit(1);
    }
    // 先寫 BMP Header (54 bytes)
    fwrite(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f);

    // 再寫像素資料 (含 padding)
    for(int i=0; i < p_bmp->Vpixels; i++){
        fwrite(p_bmp->data[i], sizeof(unsigned char), p_bmp->Hbytes, f);
    }
    fclose(f);
    return 1;
}

int bmp_free(bmp *p_bmp) 
{
    for(int i=0; i < p_bmp->Vpixels; i++){
        free(p_bmp->data[i]);
        free(p_bmp->R[i]);
        free(p_bmp->G[i]);
        free(p_bmp->B[i]);
    }
    free(p_bmp->data);
    free(p_bmp->R);
    free(p_bmp->G);
    free(p_bmp->B);

    p_bmp->data = NULL;
    p_bmp->Hpixels = 0;
    p_bmp->Vpixels = 0;
    p_bmp->HeaderInfo[0] = '\0';
    p_bmp->Hbytes = 0;
    return 1;
}
// 量化表相關
void read_quant_table(const char *fname, int qtable[8][8])
{
    FILE *fp = fopen(fname, "r");
    if(!fp){
        fprintf(stderr, "Cannot open quant table: %s\n", fname);
        exit(1);
    }
    for(int i=0; i<8; i++){
        for(int j=0; j<8; j++){
            fscanf(fp, "%d", &qtable[i][j]);
        }
    }
    fclose(fp);
}

// IDCT
static void idct_8x8(float inBlock[8][8], float outBlock[8][8])
{
    const float c0 = 1.0f / (float)sqrt(8.0);
    const float c1 = (float)sqrt(2.0/8.0);

    for(int m = 0; m < 8; m++){
        for(int n = 0; n < 8; n++){
            double sum = 0.0;
            for(int u=0; u<8; u++){
                for(int v=0; v<8; v++){
                    double cu = (u == 0) ? c0 : c1;
                    double cv = (v == 0) ? c0 : c1;
                    double val = inBlock[u][v];
                    sum += cu * cv * val
                           * cos(((2*m+1)*u*PI)/16.0)
                           * cos(((2*n+1)*v*PI)/16.0);
                }
            }
            outBlock[m][n] = (float)sum;
        }
    }
}

void do_idct_for_channel(float **Freq, float **spatial, int h, int w)
{
    int mb = h / 8;
    int nb = w / 8;
    float inBlock[8][8], outBlock[8][8];

    for(int by = 0; by < mb; by++){
        for(int bx = 0; bx < nb; bx++){
            // 取出 8x8
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    inBlock[i][j] = Freq[by*8 + i][bx*8 + j];
                }
            }
            // IDCT
            idct_8x8(inBlock, outBlock);
            // 放回
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    spatial[by*8 + i][bx*8 + j] = outBlock[i][j];
                }
            }
        }
    }
}

// 以 BT.709 公式做 YCbCr → RGB (函式)
void ycbcr_to_rgb(
    float **Y, float **Cb, float **Cr,
    unsigned char **outR, unsigned char **outG, unsigned char **outB,
    int h, int w)
{
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            float y  = Y[i][j];
            float cb = Cb[i][j] - 128.0f; 
            float cr = Cr[i][j] - 128.0f; 

            // BT.709 反向轉換
            // R = Y + 1.5748 * Cr
            // G = Y - 0.187324 * Cb - 0.468124 * Cr
            // B = Y + 1.8556 * Cb
            float r = y + 1.5748f   * cr;
            float g = y - 0.187324f * cb - 0.468124f * cr;
            float b = y + 1.8556f   * cb;

            // clamp to [0..255]
            if(r < 0) r=0; 
            if(r>255) r=255;
            if(g < 0) g=0; 
            if(g>255) g=255;
            if(b < 0) b=0; 
            if(b>255) b=255;

            outR[i][j] = (unsigned char)(r + 0.5f);
            outG[i][j] = (unsigned char)(g + 0.5f);
            outB[i][j] = (unsigned char)(b + 0.5f);
        }
    }
}

// 載入 quantized + error 檔案
void load_quantized_short(const char *fname, short **qData, int h, int w)
{
    FILE *fp = fopen(fname, "rb");
    if(!fp){
        fprintf(stderr, "Cannot open %s\n", fname);
        exit(1);
    }
    for(int i=0; i<h; i++){
        fread(qData[i], sizeof(short), w, fp);
    }
    fclose(fp);
}

void load_quantized_error(const char *fname, float **eData, int h, int w)
{
    FILE *fp = fopen(fname, "rb");
    if(!fp){
        fprintf(stderr, "Cannot open %s\n", fname);
        exit(1);
    }
    for(int i=0; i<h; i++){
        fread(eData[i], sizeof(float), w, fp);
    }
    fclose(fp);
}
int getPaddedSize(int x) {
    // 回傳最小的 >= x 的 8 倍數
    // (x + 7)/8 * 8  = ceil(x/8.0)*8
    return ((x + 7) / 8) * 8;
}
void compute_sqnr(
    unsigned char **origR, unsigned char **origG, unsigned char **origB,
    unsigned char **recR,  unsigned char **recG,  unsigned char **recB,
    int h, int w)
{
    double sumOrigR = 0.0, sumDiffR = 0.0;
    double sumOrigG = 0.0, sumDiffG = 0.0;
    double sumOrigB = 0.0, sumDiffB = 0.0;

    for(int i=0; i<h; i++){
        for(int j=0; j<w; j++){
            double oR = (double)origR[i][j];
            double rR = (double)recR[i][j];
            sumOrigR += (oR * oR);
            sumDiffR += (oR - rR)*(oR - rR);

            double oG = (double)origG[i][j];
            double rG = (double)recG[i][j];
            sumOrigG += (oG * oG);
            sumDiffG += (oG - rG)*(oG - rG);

            double oB = (double)origB[i][j];
            double rB = (double)recB[i][j];
            sumOrigB += (oB * oB);
            sumDiffB += (oB - rB)*(oB - rB);
        }
    }

    // 10*log10( sum(orig^2) / sum((orig - rec)^2) )
    double sqnrR = 10.0 * log10( sumOrigR / (sumDiffR + 1e-15) );
    double sqnrG = 10.0 * log10( sumOrigG / (sumDiffG + 1e-15) );
    double sqnrB = 10.0 * log10( sumOrigB / (sumDiffB + 1e-15) );

    printf("\nSQNR (dB):\n");
    printf("  R: %.2f dB\n", sqnrR);
    printf("  G: %.2f dB\n", sqnrG);
    printf("  B: %.2f dB\n", sqnrB);
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

// 從檔案 (ASCII 文字檔) 讀取一個 8x8 block 的 RLE，復原成 block[64] (ZigZag 順序)
void rle_decode_block_ascii(FILE *fp, short blockZ[64])
{
    // 先將 blockZ 清 0
    for(int i=0; i<64; i++){
        blockZ[i] = 0;
    }
    int count = 0;
    // 連續讀 tokens，直到填滿 64 筆為止
    while(count < 64){
        char token[16];
        if(fscanf(fp, "%s", token) != 1){
            // EOF 或 讀取失敗
            fprintf(stderr, "Error reading RLE block (ASCII)\n");
            break;
        }
        if(strncmp(token, "$skip", 5)==0){
            int num;
            fscanf(fp, "%d", &num); // 讀 skip 次數
            count += num;
        }
        else if(strncmp(token, "$value", 6)==0){
            int val;
            fscanf(fp, "%d", &val);
            blockZ[count] = (short)val;
            count++;
        }
        else {
            continue;
        }
        if(count > 64){
            fprintf(stderr, "Block overflow in ASCII RLE decode!\n");
            break;
        }
    }
}


// 讀取 ASCII RLE -> 還原成三個 channel (qFY, qFCb, qFCr)
void decode_rle_ascii(FILE *fp, int width, int height, 
                      short **qFY, short **qFCb, short **qFCr,
                      int padH, int padW)
{
    // 依 encoder 順序：逐個 block 的 Y -> Cb -> Cr
    // 但檔案裡每個 block 有三行 "(m,n, Y)  $skip.. $value.."
    //                   "(m,n, Cb) $skip.. $value.."
    //                   "(m,n, Cr) $skip.. $value.."
    int mb = padH / 8;
    int nb = padW / 8;

    // 略過檔案中的第一行 (width height)，因為在 main() 早已讀過
    // or 如果 main() 還沒讀，那此函式裡也能 fscanf，但要小心指針位置

    // ZigZag array => blockZ
    short blockZ[64];
    // 反 ZigZag => block_in
    short block_in[64];

    // DPCM reconstruction 需要記錄前一個 DC
    // => 每個 channel 分別記錄
    short prevDC_Y  = 0;
    short prevDC_Cb = 0;
    short prevDC_Cr = 0;

    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            // =========== Y ============
            // 先讀 "(by,bx,Y)" + RLE
            rle_decode_block_ascii(fp, blockZ);
            // DC reconstruction
            blockZ[0] = blockZ[0] + prevDC_Y;
            prevDC_Y  = blockZ[0];
            // 反 ZigZag => block_in
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            // 放入 qFY
            int idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFY[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }

            // =========== Cb ============
            rle_decode_block_ascii(fp, blockZ);
            blockZ[0] = blockZ[0] + prevDC_Cb;
            prevDC_Cb = blockZ[0];
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFCb[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }

            // =========== Cr ============
            rle_decode_block_ascii(fp, blockZ);
            blockZ[0] = blockZ[0] + prevDC_Cr;
            prevDC_Cr = blockZ[0];
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFCr[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }
        }
    }
}

// 讀取一個 block (64 筆) 的 RLE (binary)，復原到 blockZ[64]
void rle_decode_block_binary(FILE *fp, short blockZ[64])
{
    for(int i=0; i<64; i++){
        blockZ[i] = 0;
    }
    int count = 0;
    while(count < 64){
        unsigned char op;
        if(fread(&op, 1, 1, fp) != 1){
            fprintf(stderr, "Binary RLE read fail!\n");
            break;
        }
        short val;
        if(fread(&val, 2, 1, fp) != 1){
            fprintf(stderr, "Binary RLE read fail (val)!\n");
            break;
        }
        if(op == 0){
            // skip
            count += val;
        } else {
            // value
            blockZ[count] = val;
            count++;
        }
        if(count > 64){
            fprintf(stderr, "Block overflow in binary RLE decode!\n");
            break;
        }
    }
}
void decode_rle_binary(FILE *fp, int width, int height,
                       short **qFY, short **qFCb, short **qFCr,
                       int padH, int padW)
{
    int mb = padH / 8;
    int nb = padW / 8;

    // ZigZag
    short blockZ[64], block_in[64];

    // DPCM
    short prevDC_Y  = 0;
    short prevDC_Cb = 0;
    short prevDC_Cr = 0;

    // 先讀 Y-block 全部
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            rle_decode_block_binary(fp, blockZ);
            blockZ[0] = blockZ[0] + prevDC_Y;
            prevDC_Y = blockZ[0];
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            int idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFY[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }
        }
    }
    // 再讀 Cb-block 全部
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            rle_decode_block_binary(fp, blockZ);
            blockZ[0] = blockZ[0] + prevDC_Cb;
            prevDC_Cb = blockZ[0];
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            int idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFCb[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }
        }
    }
    // 再讀 Cr-block
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            rle_decode_block_binary(fp, blockZ);
            blockZ[0] = blockZ[0] + prevDC_Cr;
            prevDC_Cr = blockZ[0];
            for(int i=0; i<64; i++){
                block_in[ ZigZagOrder[i] ] = blockZ[i];
            }
            int idx=0;
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    qFCr[by*8 + i][bx*8 + j] = block_in[idx++];
                }
            }
        }
    }
}


int main(int argc, char **argv)
{

    int mode = atoi(argv[1]);
    // decoder 0 <outputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>
    if(mode == 0)
    {
        char *fn_output_bmp = argv[2]; 
        char *fn_r = argv[3];         
        char *fn_g = argv[4];         
        char *fn_b = argv[5];         
        char *fn_dim = argv[6];       

        bmp myBmp;
        memset(&myBmp, 0, sizeof(bmp));

        // (1) 從 dim.txt 讀取寬高
        FILE *fdim = fopen(fn_dim, "r");
        if(!fdim){
            fprintf(stderr, "Cannot open %s\n", fn_dim);
            return 1;
        }
        int width, height;
        fscanf(fdim, "%d %d", &width, &height);
        fclose(fdim);

        myBmp.Hpixels = width;
        myBmp.Vpixels = height;
        myBmp.Hbytes  = (width * 3 + 3) & (~3);

        init_bmp_header(&myBmp);

        // (2) 分配記憶體給 data, R, G, B
        myBmp.data = alloc2D_uc(height, myBmp.Hbytes);
        myBmp.R    = alloc2D_uc(height, width);
        myBmp.G    = alloc2D_uc(height, width);
        myBmp.B    = alloc2D_uc(height, width);

        // (3) 讀取 R, G, B (文字檔)
        FILE *fr = fopen(fn_r, "r");
        if(!fr){
            fprintf(stderr, "Cannot open %s\n", fn_r);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                int val;
                fscanf(fr, "%d", &val);
                myBmp.R[i][j] = (unsigned char)val;
            }
        }
        fclose(fr);

        FILE *fg = fopen(fn_g, "r");
        if(!fg){
            fprintf(stderr, "Cannot open %s\n", fn_g);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                int val;
                fscanf(fg, "%d", &val);
                myBmp.G[i][j] = (unsigned char)val;
            }
        }
        fclose(fg);

        FILE *fb = fopen(fn_b, "r");
        if(!fb){
            fprintf(stderr, "Cannot open %s\n", fn_b);
            return 1;
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                int val;
                fscanf(fb, "%d", &val);
                myBmp.B[i][j] = (unsigned char)val;
            }
        }
        fclose(fb);

        // (4) 組裝 data (BMP 預設 BGR 排列)
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                myBmp.data[i][3*j + 0] = myBmp.B[i][j]; 
                myBmp.data[i][3*j + 1] = myBmp.G[i][j];
                myBmp.data[i][3*j + 2] = myBmp.R[i][j];
            }
            // 若有 padding
            for(int p=width*3; p<myBmp.Hbytes; p++){
                myBmp.data[i][p] = 0;
            }
        }

        // (5) 寫出 BMP
        bmp_save_fn(fn_output_bmp, &myBmp);

        // (6) 釋放
        bmp_free(&myBmp);

        printf("decoder (mode 0) done.\n");
    }
    // decoder 1 <outputBMP> <Qt_Y.txt> <Qt_Cb.txt> <Qt_Cr.txt> <dim.txt>
    //            <qF_Y.raw> <qF_Cb.raw> <qF_Cr.raw> <eF_Y.raw> <eF_Cb.raw> <eF_Cr.raw>
    else if(mode == 1){
        if(argc == 13){
        char *fn_output_bmp = argv[2];
        char *fn_qtY        = argv[3];
        char *fn_qtCb       = argv[4];
        char *fn_qtCr       = argv[5];
        char *fn_dim        = argv[6];
        char *fn_qFY        = argv[7];
        char *fn_qFCb       = argv[8];
        char *fn_qFCr       = argv[9];
        char *fn_eFY        = argv[10];
        char *fn_eFCb       = argv[11];
        char *fn_eFCr       = argv[12];

        // (1) 讀取維度
        FILE *fdim = fopen(fn_dim, "r");
        if(!fdim){
            fprintf(stderr, "Cannot open %s\n", fn_dim);
            return 1;
        }
        int width, height;
        int padW, padH;
        fscanf(fdim, "%d %d", &width, &height);
        fclose(fdim);
        padW = getPaddedSize(width);
        padH = getPaddedSize(height);

        // (2) 讀取量化表
        int qtableY[8][8], qtableCb[8][8], qtableCr[8][8];
        read_quant_table(fn_qtY,  qtableY);
        read_quant_table(fn_qtCb, qtableCb);
        read_quant_table(fn_qtCr, qtableCr);

        // (3) 配置記憶體給 qF, eF
        short **qFY   = alloc2D_short(padH, padW);
        short **qFCb  = alloc2D_short(padH, padW);
        short **qFCr  = alloc2D_short(padH, padW);

        float **eFY   = alloc2D_float(padH, padW);
        float **eFCb  = alloc2D_float(padH, padW);
        float **eFCr  = alloc2D_float(padH, padW);

        // (4) 從檔案讀 qF_*.raw, eF_*.raw
        load_quantized_short(fn_qFY,  qFY,   padH, padW);
        load_quantized_short(fn_qFCb, qFCb,  padH, padW);
        load_quantized_short(fn_qFCr, qFCr,  padH, padW);

        load_quantized_error(fn_eFY,  eFY,   padH, padW);
        load_quantized_error(fn_eFCb, eFCb,  padH, padW);
        load_quantized_error(fn_eFCr, eFCr,  padH, padW);

        // (5) 重建頻率值 F = qF*Qt + eF
        float **FY   = alloc2D_float(padH, padW);
        float **FCb  = alloc2D_float(padH, padW);
        float **FCr  = alloc2D_float(padH, padW);

        int mb = padH / 8;
        int nb = padW / 8;
        for(int by=0; by<mb; by++){
            for(int bx=0; bx<nb; bx++){
                for(int i=0; i<8; i++){
                    for(int j=0; j<8; j++){
                        int row = by*8 + i;
                        int col = bx*8 + j;

                        int Qy  = qtableY[i][j];
                        int Qcb = qtableCb[i][j];
                        int Qcr = qtableCr[i][j];

                        // Y
                        float fqY = (float)qFY[row][col] * (float)Qy;
                        FY[row][col] = fqY + eFY[row][col];

                        // Cb
                        float fqCb = (float)qFCb[row][col] * (float)Qcb;
                        FCb[row][col] = fqCb + eFCb[row][col];

                        // Cr
                        float fqCr = (float)qFCr[row][col] * (float)Qcr;
                        FCr[row][col] = fqCr + eFCr[row][col];
                    }
                }
            }
        }

        // (6) 配置空間域 Y, Cb, Cr
        float **spatialY  = alloc2D_float(padH, padW);
        float **spatialCb = alloc2D_float(padH, padW);
        float **spatialCr = alloc2D_float(padH, padW);

        // (7) 做 IDCT
        do_idct_for_channel(FY,  spatialY,  padH, padW);
        do_idct_for_channel(FCb, spatialCb, padH, padW);
        do_idct_for_channel(FCr, spatialCr, padH, padW);

        // (8) 只需 [0..height-1, 0..width-1] 這區域 
        unsigned char **outR = alloc2D_uc(height, width);
        unsigned char **outG = alloc2D_uc(height, width);
        unsigned char **outB = alloc2D_uc(height, width);
        ycbcr_to_rgb(
            spatialY, 
            spatialCb, 
            spatialCr, 
            outR, 
            outG, 
            outB, 
            height, 
            width);
        // (10) 組裝 BMP
        bmp myBmp;
        memset(&myBmp, 0, sizeof(myBmp));
        myBmp.Hpixels = width;
        myBmp.Vpixels = height;
        myBmp.Hbytes  = (width*3 + 3) & (~3);

        init_bmp_header(&myBmp);

        myBmp.data = alloc2D_uc(height, myBmp.Hbytes);
        myBmp.R    = alloc2D_uc(height, width);
        myBmp.G    = alloc2D_uc(height, width);
        myBmp.B    = alloc2D_uc(height, width);
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                myBmp.R[i][j] = outR[i][j];
                myBmp.G[i][j] = outG[i][j];
                myBmp.B[i][j] = outB[i][j];
            }
        }
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                myBmp.data[i][3*j+0] = myBmp.B[i][j];
                myBmp.data[i][3*j+1] = myBmp.G[i][j];
                myBmp.data[i][3*j+2] = myBmp.R[i][j];
            }
            // padding
            for(int p=width*3; p<myBmp.Hbytes; p++){
                myBmp.data[i][p] = 0;
            }
        }

        // (11) 寫出 BMP
        bmp_save_fn(fn_output_bmp, &myBmp);

        // (12) 釋放
        free2D_short(qFY,   padH);
        free2D_short(qFCb,  padH);
        free2D_short(qFCr,  padH);

        free2D_float(eFY,   padH);
        free2D_float(eFCb,  padH);
        free2D_float(eFCr,  padH);

        free2D_float(FY,    padH);
        free2D_float(FCb,   padH);
        free2D_float(FCr,   padH);

        free2D_float(spatialY,  padH);
        free2D_float(spatialCb, padH);
        free2D_float(spatialCr, padH);

        free2D_uc(outR, height);
        free2D_uc(outG, height);
        free2D_uc(outB, height);

        bmp_free(&myBmp);

        printf("decoder (mode 1b) done. Reconstructed BMP: %s\n", fn_output_bmp);
        }
        else if(argc == 11){   //decoder 1 <QResKimberly.bmp> <Kimberly.bmp> <Qt_Y.txt> <Qt_Cb.txt> <Qt_Cr.txt> <dim.txt> <qF_Y.raw> <qF_Cb.raw> <qF_Cr.raw>
            char *fn_output_bmp = argv[2];
            char *fn_org_bmp = argv[3];
            char *fn_qtY = argv[4];
            char *fn_qtCb = argv[5];
            char *fn_qtCr = argv[6];
            char *fn_dim = argv[7];
            char *fn_qFY = argv[8];
            char *fn_qFCb = argv[9];
            char *fn_qFCr = argv[10];

            FILE *fdim = fopen(fn_dim, "r");
            if(!fdim){
                fprintf(stderr, "Cannot open %s\n", fn_dim);
                return 1;
            }
            int width, height;
            fscanf(fdim, "%d %d", &width, &height);
            fclose(fdim);

            int padW = getPaddedSize(width);
            int padH = getPaddedSize(height);

            int qtableY[8][8], qtableCb[8][8], qtableCr[8][8];
            read_quant_table(fn_qtY,  qtableY);
            read_quant_table(fn_qtCb, qtableCb);
            read_quant_table(fn_qtCr, qtableCr);

            // (3) 分配記憶體給 qF
            short **qFY   = alloc2D_short(padH, padW);
            short **qFCb  = alloc2D_short(padH, padW);
            short **qFCr  = alloc2D_short(padH, padW);

            // (4) load qF_*.raw
            load_quantized_short(fn_qFY,  qFY,   padH, padW);
            load_quantized_short(fn_qFCb, qFCb,  padH, padW);
            load_quantized_short(fn_qFCr, qFCr,  padH, padW);

            // (5) 重建頻率值 F = qF * Qt (無誤差檔)
            float **FY   = alloc2D_float(padH, padW);
            float **FCb  = alloc2D_float(padH, padW);
            float **FCr  = alloc2D_float(padH, padW);

            int mb = padH / 8;
            int nb = padW / 8;
            for(int by=0; by<mb; by++){
                for(int bx=0; bx<nb; bx++){
                    for(int i=0; i<8; i++){
                        for(int j=0; j<8; j++){
                            int row = by*8 + i;
                            int col = bx*8 + j;
                            FY[row][col]   = (float)qFY[row][col]  * (float)qtableY[i][j];
                            FCb[row][col]  = (float)qFCb[row][col] * (float)qtableCb[i][j];
                            FCr[row][col]  = (float)qFCr[row][col] * (float)qtableCr[i][j];
                        }
                    }
                }
            }

            float **spatialY  = alloc2D_float(padH, padW);
            float **spatialCb = alloc2D_float(padH, padW);
            float **spatialCr = alloc2D_float(padH, padW);
            do_idct_for_channel(FY,  spatialY,  padH, padW);
            do_idct_for_channel(FCb, spatialCb, padH, padW);
            do_idct_for_channel(FCr, spatialCr, padH, padW);

            // (7) 只需 [0..height-1, 0..width-1]
            unsigned char **outR = alloc2D_uc(height, width);
            unsigned char **outG = alloc2D_uc(height, width);
            unsigned char **outB = alloc2D_uc(height, width);

            ycbcr_to_rgb(spatialY, spatialCb, spatialCr, outR, outG, outB, height, width);

            // (8) 組裝 BMP
            bmp myBmp;
            memset(&myBmp, 0, sizeof(myBmp));
            myBmp.Hpixels = width;
            myBmp.Vpixels = height;
            myBmp.Hbytes  = (width*3 + 3) & (~3);
            init_bmp_header(&myBmp);

            myBmp.data = alloc2D_uc(height, myBmp.Hbytes);
            myBmp.R    = alloc2D_uc(height, width);
            myBmp.G    = alloc2D_uc(height, width);
            myBmp.B    = alloc2D_uc(height, width);

            for(int i=0; i<height; i++){
                for(int j=0; j<width; j++){
                    myBmp.R[i][j] = outR[i][j];
                    myBmp.G[i][j] = outG[i][j];
                    myBmp.B[i][j] = outB[i][j];
                }
            }
            for(int i=0; i<height; i++){
                for(int j=0; j<width; j++){
                    myBmp.data[i][3*j+0] = myBmp.B[i][j];
                    myBmp.data[i][3*j+1] = myBmp.G[i][j];
                    myBmp.data[i][3*j+2] = myBmp.R[i][j];
                }
                for(int p=width*3; p<myBmp.Hbytes; p++){
                    myBmp.data[i][p] = 0;
                }
            }
            // (9) 寫出重建後 BMP
            bmp_save_fn(fn_output_bmp, &myBmp);

            // (10) 讀取原圖 BMP，並計算 SQNR
            bmp origBmp;
            memset(&origBmp, 0, sizeof(origBmp));
            if(!bmp_load_fn(fn_org_bmp, &origBmp)){
                fprintf(stderr, "Load original BMP fail!\n");
            } else {
                // 確定大小一致
                if(origBmp.Hpixels == width && origBmp.Vpixels == height){
                    compute_sqnr(origBmp.R, origBmp.G, origBmp.B, 
                                 myBmp.R,  myBmp.G,  myBmp.B, 
                                 height, width);
                } else {
                    fprintf(stderr, "Warning: dimension mismatch. Can't compute SQNR!\n");
                }
                bmp_free(&origBmp);
            }

            // (11) 釋放
            free2D_short(qFY,   padH);
            free2D_short(qFCb,  padH);
            free2D_short(qFCr,  padH);
            free2D_float(FY,    padH);
            free2D_float(FCb,   padH);
            free2D_float(FCr,   padH);
            free2D_float(spatialY,  padH);
            free2D_float(spatialCb, padH);
            free2D_float(spatialCr, padH);

            free2D_uc(outR, height);
            free2D_uc(outG, height);
            free2D_uc(outB, height);

            bmp_free(&myBmp);

            printf("decoder (mode 1a) done.");
        }
        else
        {
            fprintf(stderr, "Not enough (or too many) arguments for mode=1\n");
            return 1;
        }
    }
    
    else if(mode == 2)
{
    // decoder 2 <outputBMP> <ascii|binary> <rle_file>
    //    => 從 rle_file 讀取 RLE, 還原 Y/Cb/Cr => IDCT => BMP
    if(argc < 5){
        fprintf(stderr, "Usage: decoder 2 <outputBMP> <ascii|binary> <rle_file>\n");
        return 0;
    }
    char *fn_output_bmp = argv[2];
    char *fmt_str       = argv[3]; // "ascii" or "binary"
    char *fn_rle        = argv[4];

    // 1) 打開 RLE 檔
    FILE *fp_rle = NULL;
    if(strcmp(fmt_str,"binary")==0){
        fp_rle = fopen(fn_rle, "rb");
    } else {
        fp_rle = fopen(fn_rle, "r");
    }
    if(!fp_rle){
        fprintf(stderr, "Cannot open RLE file: %s\n", fn_rle);
        return 1;
    }

    // 2) 先讀取圖像大小 (width, height)
    int width, height;
    if(strcmp(fmt_str,"ascii")==0){
        // ASCII 版: 第一行是 "width height"
        fscanf(fp_rle, "%d %d", &width, &height);
    } else {
        // Binary 版: 開頭 8 bytes = (int)width, (int)height
        fread(&width,  sizeof(int), 1, fp_rle);
        fread(&height, sizeof(int), 1, fp_rle);
    }

    // 3) 計算 padding
    int padW = getPaddedSize(width);
    int padH = getPaddedSize(height);

    // 4) 分配 qF (2D short) 給三通道
    short **qFY  = alloc2D_short(padH, padW);
    short **qFCb = alloc2D_short(padH, padW);
    short **qFCr = alloc2D_short(padH, padW);

    // 5) 依 ascii/binary 讀 RLE => 填入 qFY, qFCb, qFCr
    if(strcmp(fmt_str,"ascii")==0){
        decode_rle_ascii(fp_rle, width, height, qFY, qFCb, qFCr, padH, padW);
    } else {
        decode_rle_binary(fp_rle, width, height, qFY, qFCb, qFCr, padH, padW);
    }
    fclose(fp_rle);

    // 6) 反量化 => F[u,v] = qF[u,v] * Qtable
    //    (跟 encoder 裡的 quant 相反)
    float **FY  = alloc2D_float(padH, padW);
    float **FCb = alloc2D_float(padH, padW);
    float **FCr = alloc2D_float(padH, padW);

    int mb = padH/8, nb = padW/8;
    for(int by=0; by<mb; by++){
        for(int bx=0; bx<nb; bx++){
            for(int i=0; i<8; i++){
                for(int j=0; j<8; j++){
                    int row = by*8 + i;
                    int col = bx*8 + j;
                    FY[row][col]  = (float)qFY[row][col]  * (float)qtableY[i][j];
                    FCb[row][col] = (float)qFCb[row][col] * (float)qtableCb[i][j];
                    FCr[row][col] = (float)qFCr[row][col] * (float)qtableCr[i][j];
                }
            }
        }
    }

    // 7) IDCT => 得到空間域 Y, Cb, Cr
    float **spatialY  = alloc2D_float(padH, padW);
    float **spatialCb = alloc2D_float(padH, padW);
    float **spatialCr = alloc2D_float(padH, padW);

    do_idct_for_channel(FY,  spatialY,  padH, padW);
    do_idct_for_channel(FCb, spatialCb, padH, padW);
    do_idct_for_channel(FCr, spatialCr, padH, padW);

    // 8) 製作 RGB
    unsigned char **outR = alloc2D_uc(height, width);
    unsigned char **outG = alloc2D_uc(height, width);
    unsigned char **outB = alloc2D_uc(height, width);

    // 只取 [0..height-1, 0..width-1]
    ycbcr_to_rgb(spatialY, spatialCb, spatialCr, outR, outG, outB, height, width);

    // 9) 組裝 BMP
    bmp myBmp;
    memset(&myBmp, 0, sizeof(myBmp));
    myBmp.Hpixels = width;
    myBmp.Vpixels = height;
    myBmp.Hbytes  = (width*3 + 3) & (~3);
    init_bmp_header(&myBmp);

    myBmp.data = alloc2D_uc(height, myBmp.Hbytes);
    myBmp.R    = alloc2D_uc(height, width);
    myBmp.G    = alloc2D_uc(height, width);
    myBmp.B    = alloc2D_uc(height, width);

    // 塞回 R,G,B
    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){
            myBmp.R[i][j] = outR[i][j];
            myBmp.G[i][j] = outG[i][j];
            myBmp.B[i][j] = outB[i][j];
        }
    }
    // data = BGR
    for(int i=0; i<height; i++){
        for(int j=0; j<width; j++){
            myBmp.data[i][3*j + 0] = myBmp.B[i][j];
            myBmp.data[i][3*j + 1] = myBmp.G[i][j];
            myBmp.data[i][3*j + 2] = myBmp.R[i][j];
        }
        // 填充 padding
        for(int p=width*3; p<myBmp.Hbytes; p++){
            myBmp.data[i][p] = 0;
        }
    }

    // 10) 寫出 BMP
    bmp_save_fn(fn_output_bmp, &myBmp);

    // 11) 釋放
    free2D_short(qFY,  padH);
    free2D_short(qFCb, padH);
    free2D_short(qFCr, padH);
    free2D_float(FY,  padH);
    free2D_float(FCb, padH);
    free2D_float(FCr, padH);
    free2D_float(spatialY,  padH);
    free2D_float(spatialCb, padH);
    free2D_float(spatialCr, padH);
    free2D_uc(outR, height);
    free2D_uc(outG, height);
    free2D_uc(outB, height);
    bmp_free(&myBmp);

    printf("decoder (mode 2) done.");
}
    
    else {
        fprintf(stderr, "Unknown mode!\n");
    }
    return 0;
}
