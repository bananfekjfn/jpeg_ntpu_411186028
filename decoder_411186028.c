//   decoder 0 <outputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>
//  1. 從 dim.txt 取得寬高
//  2. 讀取 R, G, B 三個文字檔
//  3. 用「原封不動的 Header」 + (新的) 像素資料，產生 outputBMP
//  => 只要 R/G/B 未被改動，就能 bit-for-bit 重現原 BMP
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _bmp
{
    int Hpixels; // 圖片寬度
    int Vpixels; 
    unsigned char HeaderInfo[54];
    unsigned long int Hbytes;
    unsigned char **data; // 二維陣列儲存像素資料
    unsigned char **R;  //二維陣列儲存 R 
    unsigned char **G;
    unsigned char **B;
    float **Y;    
    float **Cb;
    float **Cr;
} bmp;

void init_bmp_header(bmp *p_bmp) {  // 初始化 BMP 的 Header
    memset(p_bmp->HeaderInfo, 0, 54);
    p_bmp->HeaderInfo[0] = 'B';
    p_bmp->HeaderInfo[1] = 'M';
    // 設定 BMP 檔案大小
    int fileSize = (p_bmp->Hbytes * p_bmp->Vpixels) + 54;
    memcpy(&p_bmp->HeaderInfo[2], &fileSize, 4);

    int bfOffBits = 54;// 資料的起始偏移量
    memcpy(&p_bmp->HeaderInfo[10], &bfOffBits, 4);

    int biSize = 40;
    memcpy(&p_bmp->HeaderInfo[14], &biSize, 4);

    memcpy(&p_bmp->HeaderInfo[18], &p_bmp->Hpixels, 4);
    memcpy(&p_bmp->HeaderInfo[22], &p_bmp->Vpixels, 4);

    short biPlanes = 1; // 平面數，固定為 1
    memcpy(&p_bmp->HeaderInfo[26], &biPlanes, 2);

    short biBitCount = 24;   // 每個像素的位數
    memcpy(&p_bmp->HeaderInfo[28], &biBitCount, 2);

    int biSizeImage = p_bmp->Hbytes * p_bmp->Vpixels;
    memcpy(&p_bmp->HeaderInfo[34], &biSizeImage, 4);
    // X/Y 方向的像素密度
    int biXPelsPerMeter = 2835; 
    memcpy(&p_bmp->HeaderInfo[38], &biXPelsPerMeter, 4);

    int biYPelsPerMeter = 2835;
    memcpy(&p_bmp->HeaderInfo[42], &biYPelsPerMeter, 4);
}

// 直接用原來的 HeaderInfo[0..53]，再把像素資料寫到檔案。
// 若同時保有相同 R/G/B，則生成之 BMP 應與原始檔 bit-for-bit 相同。
int bmp_save_fn(char *filename, bmp *p_bmp)  // 儲存 BMP 檔案
{
    FILE* f = fopen(filename, "wb");
    if(f == NULL){
        printf("\n\nFILE CREATION ERROR: %s\n\n", filename);
        exit(1);
    }
    // 先寫入原封不動的 BMP Header(54 bytes)
    fwrite(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f);

    // 再寫像素資料 (含 padding)
    for(int i=0; i < p_bmp->Vpixels; i++){
        fwrite(p_bmp->data[i], sizeof(unsigned char), p_bmp->Hbytes, f);
    }

    fclose(f);
    return 1;
}

// 釋放 bmp 的動態配置記憶體
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

int main(int argc, char **argv)
{
    if(argc < 7){
        fprintf(stderr, "\nUsage (mode 0):\n");
        fprintf(stderr, "  decoder 0 <outputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>\n");
        return 0;
    }
//   decoder 0 <outputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>
    int mode = atoi(argv[1]);
    char *fn_output_bmp = argv[2]; // 產生的 BMP，e.g. ResKimberly.bmp
    char *fn_r = argv[3];         // R.txt
    char *fn_g = argv[4];         //  G.txt
    char *fn_b = argv[5];         //  B.txt
    char *fn_dim = argv[6];       //  dim.txt

    if(mode == 0)
    {
        bmp myBmp;
        memset(&myBmp, 0, sizeof(bmp));

        // 1) 從 dim.txt 讀取寬高 
        FILE *fdim = fopen(fn_dim, "r");
        if(!fdim){
            fprintf(stderr, "Cannot open %s\n", fn_dim);
            return 1;
        }

        int width, height;
        fscanf(fdim, "%d %d", &width, &height);  // 第一行：width, height
        myBmp.Hpixels = width;
        myBmp.Vpixels = height;
        fclose(fdim);

        // 通常 (width * 3 + 3) & (~3) 和原 Header 應該相同，否則 BMP 也不會 bit for bit 相同
        int RowBytes = (width * 3 + 3) & (~3);
        myBmp.Hbytes = RowBytes;

        init_bmp_header(&myBmp);

        // 2) 分配記憶體給 data, R, G, B
        myBmp.data = (unsigned char **)malloc(height * sizeof(unsigned char*));
        myBmp.R = (unsigned char **)malloc(height * sizeof(unsigned char*));
        myBmp.G = (unsigned char **)malloc(height * sizeof(unsigned char*));
        myBmp.B = (unsigned char **)malloc(height * sizeof(unsigned char*));

        for(int i=0; i<height; i++){
            myBmp.data[i] = (unsigned char*)malloc(RowBytes * sizeof(unsigned char));
            myBmp.R[i]    = (unsigned char*)malloc(width * sizeof(unsigned char));
            myBmp.G[i]    = (unsigned char*)malloc(width * sizeof(unsigned char));
            myBmp.B[i]    = (unsigned char*)malloc(width * sizeof(unsigned char));
        }

        // 3) 讀取 R, G, B 文字檔，覆寫 myBmp.R/G/B
        // --- R ---
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

        // --- G ---
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

        // --- B ---
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

        // 4) 將 R, G, B 三個 channel 寫回去來組裝 data (BMP預設 BGR)
        for(int i=0; i<height; i++){
            for(int j=0; j<width; j++){
                myBmp.data[i][3*j + 0] = myBmp.B[i][j]; // B
                myBmp.data[i][3*j + 1] = myBmp.G[i][j]; // G
                myBmp.data[i][3*j + 2] = myBmp.R[i][j]; // R
            }
            // 如果 RowBytes > width*3，就要補 0 做 padding
            for(int p = width*3; p < RowBytes; p++){
                myBmp.data[i][p] = 0;
            }
        }

        // 5) 寫出 BMP：HeaderInfo (原封不動) + data
        bmp_save_fn(fn_output_bmp, &myBmp);

        // 6) 釋放記憶體
        bmp_free(&myBmp);

        printf("decoder (mode 0) done.\n");
    }
    else {
        fprintf(stderr, "Only mode 0 is implemented.\n");
    }

    return 0;
}
