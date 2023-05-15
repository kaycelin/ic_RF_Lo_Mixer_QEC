# 0302_RF_Lo_Mixer_QEC

## Summary of simulation
<img src="https://user-images.githubusercontent.com/87049112/226227322-193bd0cb-974c-423a-a995-c5f428c9c972.png" width="50%">
- Summary_NR Test model.md


## Main: NRTestModelWaveformGeneration_main_k
  - NR-TM or PDSCH FRC waveform parameters
  ```js
  flag_FreqRange = 'FR1' 
  dlnrref = "NR-FR1-TM3.1";   
  bw = "100MHz";
  scs = "30kHz";
  dm = "FDD";
  ncellid = 1;  % NCellID
  sv      = "16.7.0";  % TS 38.141-x version (NR-TM only)
  Nframes = 1
  ```
