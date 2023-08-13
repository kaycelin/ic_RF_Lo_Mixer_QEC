## 2023-08-12, Example of evm comparsion between TM1.1(QPSK) and TM3.1(64QAM) 
**Waveform: 100MHz scs:30kHz sampleRate:983.04MHz**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/e2fa3b33-c53a-464b-af4c-7a358811c378)

**CCDF**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/da837570-4305-4c38-bb79-64794d82d689)

**System link: Upconversion/Lo phase noise/IQ imbalance/IQ offset/Lo leakage/Downconversion/Demodulation**
- lo parameters
```js
LoParms = 

  struct with fields:

     flag_LoPhaseNoise: 1
     phsNzFreqOffsetHz: [1×9 double]
           phsNzLvldBc: [-90 -95 -100 -105 -110 -115 -120 -125 -130]
                    fs: 983040000
                   fLo: 245760000
                Nsamps: 983040
    flag_LoQuad_method: 'LoIQ->LoPN'
                PwrdBm: 0
                  fnum: 41802
           flag_LoQuad: 'IQ'

 LoParms.phsNzFreqOffsetHz = [100000, 200000, 400000, 600000, 800000, 1200000, 1800000, 6000000, 10000000]
```
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/c98989e3-5443-4850-af51-3b22b79ccb52)
         
- mixer parameters
```js
mixParms = 

  struct with fields:

                            x: [983040×1 double]
                           lo: [983040×2 double]
                           fs: 983040000
                         fnum: 41901
               flag_UDconvert: 'U'
    flag_UDconvert_FreqSelect: 'H'
                 flag_MixType: 'IQ'
             flag_IQimbalance: 1
                    imb_AmpDB: -0.5
                   imb_PhsDeg: 1
             flag_IQLvlOffset: 1
              imb_LvlOffsetDB: 0.5
                     flag_Amp: 1
              ampIn_op1dB_dBm: 35
               ampIn_oip3_dBm: 5
                ampIn_gain_dB: 10
                 ampOut_nf_dB: 5
                  ampOut_flat: 0
                    flag_Amp2: 0
                mixerOut_flat: []
               flag_LoLeakage: 1
                    Lo2IF_dBc: 60
                    Lo2RF_dBc: 60
                         fNCO: 0
              flag_IQimb_comp: 0
               imb_AmpDB_comp: 0
              imb_PhsDeg_comp: 0
        flag_IQLvlOffset_comp: 0
```

**Waveform Upconversion**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/6b86ca9e-d931-4960-a48e-cf5dcc1266bd)


**Receiver Waveform Demodulation EVM Result**
- TM3.1      
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/ca44ed66-f1e0-4932-9552-78059dd6f2a3)

- TM1.1        
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/755be181-d4ee-406b-a273-12c9eebabd57)


**Receiver Waveform Demodulation Constellation Result**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/fbe77a5f-0488-44ac-8d99-780ee4ae33ab)

**Receiver Waveform CCDF Result**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/9d22b66e-1867-4768-a089-61806bfb20e9)













