# 使用 linux 指令 diff/cmp 證明 Kimberly.bmp 與重建的圖檔 ResKimberly.bmp 是相同的！

# 編譯器與選項設定
CC = gcc
CFLAGS = -Wall -O2

# 原始碼檔案
ENCODER_SRC = encoder_411186028.c
DECODER_SRC = decoder_411186028.c

# 執行檔
ENCODER_EXE = encoder_411186028.exe
DECODER_EXE = decoder_411186028.exe

# 輸入與輸出檔案
INPUT_BMP = Kimberly.bmp
OUTPUT_BMP = ResKimberly.bmp
R_TXT = R.txt
G_TXT = G.txt
B_TXT = B.txt
DIM_TXT = dim.txt

# 預設目標：編譯所有執行檔
all: $(ENCODER_EXE) $(DECODER_EXE)

# 編譯 encoder
$(ENCODER_EXE): $(ENCODER_SRC)
	$(CC) $(CFLAGS) -o $@ $<

# 編譯 decoder
$(DECODER_EXE): $(DECODER_SRC)
	$(CC) $(CFLAGS) -o $@ $<

# demo0: 執行 mode=0 的編碼與解碼，並比較檔案
demo0: all
	@echo "Running encoder mode=0..."
	./$(ENCODER_EXE) 0 $(INPUT_BMP) $(R_TXT) $(G_TXT) $(B_TXT) $(DIM_TXT)
	@echo "Encoder completed. Running decoder..."
	./$(DECODER_EXE) 0 $(OUTPUT_BMP) $(R_TXT) $(G_TXT) $(B_TXT) $(DIM_TXT)
	@echo "Decoder completed. Comparing original and reconstructed BMP files..."
	
	# 使用 diff 比較檔案
	@diff $(INPUT_BMP) $(OUTPUT_BMP) > /dev/null && echo "Files are identical (diff)." || echo "Files are different (diff)."
	
	# 使用 cmp 比較檔案
	@cmp -s $(INPUT_BMP) $(OUTPUT_BMP) && echo "Files are identical (cmp)." || echo "Files are different (cmp)."

# demo1: 執行 mode=1 的編碼與解碼，並比較檔案
demo1: all
	@echo "Running encoder mode=1..."
	./$(ENCODER_EXE) 1 Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw
	
	@echo "Running decoder mode=1a..."
	./$(DECODER_EXE) 1 QResKimberly.bmp Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw

	@echo "Running decoder mode=1b..."
	./$(DECODER_EXE) 1 $(OUTPUT_BMP) Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw

	@echo "Comparing files using diff/cmp..."
	@diff $(INPUT_BMP) $(OUTPUT_BMP) > /dev/null && echo "Files are identical (diff)." || echo "Files are different (diff)."
	@cmp -s $(INPUT_BMP) $(OUTPUT_BMP) && echo "Files are identical (cmp)." || echo "Files are different (cmp)."
	@echo "demo1 done."

# 清除所有生成的檔案
clean:
	cmd /c "del /Q \
	$(ENCODER_EXE) \
	$(DECODER_EXE) \
	$(R_TXT) \
	$(G_TXT) \
	$(B_TXT) \
	$(DIM_TXT) \
	$(OUTPUT_BMP) \
	QResKimberly.bmp \
	ResKimberly.bmp \
	Qt_Y.txt \
	Qt_Cb.txt \
	Qt_Cr.txt \
	qF_Y.raw \
	qF_Cb.raw \
	qF_Cr.raw \
	eF_Y.raw \
	eF_Cb.raw \
	eF_Cr.raw"

	@echo "Cleaned all generated files."

# 定義 phony targets，避免與檔案名稱衝突
.PHONY: all demo0 demo1 clean
