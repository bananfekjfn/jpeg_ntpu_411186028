# jpeg_ntpu_411186028

本專案實作一套影像壓縮與解壓縮系統，模擬 JPEG 的核心流程，並以 BMP 格式圖片為輸入與輸出。整體架構涵蓋 編碼端 (encoder) 與 解碼端 (decoder)，功能包含：  
**BMP 圖片處理**   
  1.支援 BMP 格式讀取與儲存，能將原圖像資料拆解成 R/G/B 三個通道，或重建為完整的 BMP 圖片。  
  2.提供記憶體配置、釋放等通用函式，確保影像處理過程穩定。  
**色彩空間轉換 (RGB ↔ YCbCr)**  
  依據 BT.709 標準將 RGB 轉換為 Y、Cb、Cr 分量，有利於後續壓縮。  
**頻域轉換與量化**  
  1.對每個通道進行區塊化（8×8 block），並施加 2D-DCT (離散餘弦轉換)，將空間域轉換至頻率域。  
  2.使用標準 JPEG 量化表 (Y、Cb、Cr) 進行量化，並保存量化誤差。  
**熵編碼 (DPCM、ZigZag、RLE)**  
  1.對 DC 係數使用 DPCM 差分編碼。  
  2.對 8×8 頻率係數採用 ZigZag 掃描，將二維轉換為一維序列。  
  3.進一步以 RLE (Run-Length Encoding) 進行壓縮，支援 ASCII 與二進位格式。  
**解碼流程**  
  1.從壓縮檔讀回 RLE，還原 ZigZag 與 DPCM。  
  2.反量化後，使用 IDCT 將頻率域轉回空間域。  
  3.透過 YCbCr → RGB 轉換組裝回彩色圖片，並輸出為 BMP 檔案。  
**效能評估**  
  1.提供 SQNR (Signal-to-Quantization-Noise Ratio) 計算，用於量化不同頻率係數的訊號與誤差比。  
  2.計算並顯示各通道與整體的 壓縮率 (Compression Ratio)，比較原始 BMP 與壓縮檔案大小。  
