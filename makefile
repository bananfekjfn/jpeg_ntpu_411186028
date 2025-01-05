# 編譯器與編譯選項設定
CC = gcc
CFLAGS = -Wall -O2

# 原始碼檔案
ENCODER_SRC = encoder_411186028.c
DECODER_SRC = decoder_411186028.c

# 產生的執行檔
ENCODER_EXE = encoder_411186028.exe
DECODER_EXE = decoder_411186028.exe

# 輸入與輸出檔案設定
INPUT_BMP = Kimberly.bmp
OUTPUT_BMP = ResKimberly.bmp
R_TXT = R.txt
G_TXT = G.txt
B_TXT = B.txt
DIM_TXT = dim.txt

# 預設目標：編譯所有執行檔
all: $(ENCODER_EXE) $(DECODER_EXE)

# 編譯 encoder.exe
$(ENCODER_EXE): $(ENCODER_SRC)
	$(CC) $(CFLAGS) -o $@ $<

# 編譯 decoder.exe
$(DECODER_EXE): $(DECODER_SRC)
	$(CC) $(CFLAGS) -o $@ $<

# 執行 demo0：編譯、執行編碼與解碼，並比較檔案
demo0: all
	@echo Running encoder...
	./$(ENCODER_EXE) 0 $(INPUT_BMP) $(R_TXT) $(G_TXT) $(B_TXT) $(DIM_TXT)
	@echo Encoder completed. Running decoder...
	./$(DECODER_EXE) 0 $(OUTPUT_BMP) $(R_TXT) $(G_TXT) $(B_TXT) $(DIM_TXT)
	@echo Decoder completed. Comparing original and reconstructed BMP files...

	# 使用 diff 比較文件
	@diff $(INPUT_BMP) $(OUTPUT_BMP) > /dev/null && echo "Files are identical." || echo "Files are different."

	# 使用 cmp 比較文件
	@cmp -s $(INPUT_BMP) $(OUTPUT_BMP) && echo "Files are identical (cmp)." || echo "Files are different (cmp)."

# 清除所有編譯產物及生成的檔案
clean:
	cmd /c "del /Q $(ENCODER_EXE) $(DECODER_EXE) $(R_TXT) $(G_TXT) $(B_TXT) $(DIM_TXT) $(OUTPUT_BMP)"
	@echo Cleaned all generated files.

# 定義 phony 目標，避免與檔案名稱衝突
.PHONY: all demo0 clean
