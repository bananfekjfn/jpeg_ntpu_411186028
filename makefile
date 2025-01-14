CC = gcc
CFLAGS = -Wall -O2

ENCODER_SRC = encoder_411186028.c
DECODER_SRC = decoder_411186028.c
ENCODER_EXE = encoder_411186028.exe
DECODER_EXE = decoder_411186028.exe

INPUT_BMP = Kimberly.bmp
OUTPUT_BMP = ResKimberly.bmp
OUTPUTQ_BMP = QResKimberly.bmp
R_TXT = R.txt
G_TXT = G.txt
B_TXT = B.txt
DIM_TXT = dim.txt

all: $(ENCODER_EXE) $(DECODER_EXE)

$(ENCODER_EXE): $(ENCODER_SRC)
	$(CC) $(CFLAGS) -o $@ $< -lm

$(DECODER_EXE): $(DECODER_SRC)
	$(CC) $(CFLAGS) -o $@ $< -lm

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

demo1: all
	@echo "Running encoder mode=1..."
	./$(ENCODER_EXE) 1 Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw
	
	@echo "Running decoder mode=1a..."
	./$(DECODER_EXE) 1 QResKimberly1.bmp Kimberly.bmp Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw

	@echo "Running decoder mode=1b..."
	./$(DECODER_EXE) 1 $(OUTPUT_BMP) Qt_Y.txt Qt_Cb.txt Qt_Cr.txt dim.txt qF_Y.raw qF_Cb.raw qF_Cr.raw eF_Y.raw eF_Cb.raw eF_Cr.raw

	@echo "Comparing files using diff/cmp..."
	@diff $(INPUT_BMP) $(OUTPUT_BMP) > /dev/null && echo "Files are identical (diff)." || echo "Files are different (diff)."
	@cmp -s $(INPUT_BMP) $(OUTPUT_BMP) && echo "Files are identical (cmp)." || echo "Files are different (cmp)."
	@echo "demo1 done."

demo2: all
	@echo "Running encoder mode=2a..."
	./$(ENCODER_EXE) 2 Kimberly.bmp ascii rle_code.txt
	@echo ""
	@echo "Running decoder mode=2a..." 
	./$(DECODER_EXE) 2 QResKimberly2a.bmp ascii rle_code.txt 
	@echo ""
	@echo "Comparing files QResKimberly1 with QResKimberly2a using diff/cmp..."
	@cmp -s QResKimberly2a QResKimberly1 && echo "Files are identical (cmp)." || echo "Files are different (cmp)."
	@echo "demo 2a done."

	@echo "Running encoder mode=2b..."
	./$(ENCODER_EXE) 2 Kimberly.bmp binary rle_code.bin 
	@echo ""
	@echo "Running decoder mode=2b..."
	./$(DECODER_EXE) 2 QResKimberly2b.bmp binary rle_code.bin
	@echo ""
	@echo "Comparing files QResKimberly1 with QResKimberly2b using diff/cmp..."
	@cmp -s QResKimberly2b QResKimberly1 && echo "Files are identical (cmp)." || echo "Files are different (cmp)."
	@echo "demo 2b done." 

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
	rle_code.bin \
	rle_code.txt \
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
