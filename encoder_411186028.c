// Usage (mode 0):
 //  encoder 0 <inputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>
// 1. 讀入 BMP (含 54 bytes Header)
// 2. 將寬高寫到 dim.txt
// 3. 將 R, G, B channel 寫到 R.txt, G.txt, B.txt
//  => 保證未來可做 bit-for-bit 還原

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct _bmp
{
    int Hpixels; // 圖像寬度（像素數）
    int Vpixels; // 圖像高度（像素數）
    unsigned char HeaderInfo[54];   // BMP的前54 bytes
    unsigned long int Hbytes;       // 一列(含padding)的bytes數
    unsigned char **data;           // 原始 BGR 資料(含padding)
    unsigned char **R;
    unsigned char **G;
    unsigned char **B;
    float **Y;    // 以下暫不使用
    float **Cb;
    float **Cr;
} bmp;

// 讀取 BMP：讀完之後，p_bmp->HeaderInfo 會保有前 54 bytes，
// 並以 p_bmp->data[][] 儲存所有像素(含padding)，再分離 R, G, B。
int bmp_load_fn(char *filename, bmp *p_bmp)
{
    FILE* f = fopen(filename, "rb"); //以二進位讀取模式（"rb"）打開 BMP 文件，成功打開後，文件指針 f 用於後續讀取
    if(f == NULL){
        printf("\n\n%s NOT FOUND\n\n", filename);
        exit(1);
    }
    fread(p_bmp->HeaderInfo, sizeof(unsigned char), 54, f); // 讀取 BMP header (固定 54 bytes)
    
    int width, height;
    memcpy(&width, &p_bmp->HeaderInfo[18], sizeof(int));
    memcpy(&height, &p_bmp->HeaderInfo[22], sizeof(int)); // 從 BMP header 取出寬度與高度 (offset=18,22)
    
    p_bmp->Hpixels = width;
    p_bmp->Vpixels = height;

    int RowBytes = (width * 3 + 3) & (~3); //每行包含 padding 的總字節數
    p_bmp->Hbytes = RowBytes;
    
    p_bmp->data = (unsigned char **)malloc(height * sizeof(unsigned char*)); // 分配 data[][]，儲存所有的 (B, G, R, padding)
    for(int i=0; i<height; i++){
        p_bmp->data[i] = (unsigned char *)malloc(RowBytes * sizeof(unsigned char));
    }
    
    for(int i=0; i<height; i++){ // 讀取像素資料
        fread(p_bmp->data[i], sizeof(unsigned char), RowBytes, f);
    }
    fclose(f);

    // 分配並填入 R, G, B
    p_bmp->R = (unsigned char **)malloc(height * sizeof(unsigned char*));
    p_bmp->G = (unsigned char **)malloc(height * sizeof(unsigned char*));
    p_bmp->B = (unsigned char **)malloc(height * sizeof(unsigned char*));

    for(int i=0; i<height; i++){
        p_bmp->R[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        p_bmp->G[i] = (unsigned char*)malloc(width * sizeof(unsigned char));
        p_bmp->B[i] = (unsigned char*)malloc(width * sizeof(unsigned char));

        for(int j=0; j<width; j++){
            // BMP預設的像素順序為 BGR 提取 B, G, R 資料
            p_bmp->B[i][j] = p_bmp->data[i][3*j + 0];
            p_bmp->G[i][j] = p_bmp->data[i][3*j + 1];
            p_bmp->R[i][j] = p_bmp->data[i][3*j + 2];
        }
    }

    return 1; // success
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
        fprintf(stderr, "  encoder 0 <inputBMP> <R.txt> <G.txt> <B.txt> <dim.txt>\n");
        return 0;
    }
    int mode = atoi(argv[1]);
    char *fn_input_bmp = argv[2];  //  Kimberly.bmp
    char *fn_r = argv[3];         //  R.txt
    char *fn_g = argv[4];         //  G.txt
    char *fn_b = argv[5];         //  B.txt
    char *fn_dim = argv[6];       //  dim.txt
    if(mode == 0)
    {
        bmp myBmp;
        bmp_load_fn(fn_input_bmp, &myBmp);
        int width = myBmp.Hpixels;
        int height = myBmp.Vpixels;
        //--- 1) 將寬高形式寫到 dim.txt ---
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

        //3) 將 G channel 寫到 G.txt
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

        //4) 將 B channel 寫到 B.txt 
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

        //釋放記憶體 
        bmp_free(&myBmp);

        printf("encoder (mode 0) done.\n");
    }
    else {
        fprintf(stderr, "Only mode 0 is implemented.\n");
    }

    return 0;
}
