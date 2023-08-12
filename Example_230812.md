## 2023-08-12, Example of evm comparsion between TM1.1 and TM3.1 
**Waveform: 100MHz scs:30kHz sampleRate:983.04MHz**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/f97b7934-a55a-438c-beb2-f280c9583e77)
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
                    imb_AmpDB: -1
                   imb_PhsDeg: 1.5
             flag_IQLvlOffset: 1
              imb_LvlOffsetDB: 3.1
                     flag_Amp: 1
              ampIn_op1dB_dBm: 35
               ampIn_oip3_dBm: 5
                ampIn_gain_dB: 10
                 ampOut_nf_dB: 5
                  ampOut_flat: 0
                    flag_Amp2: 0
                mixerOut_flat: []
               flag_LoLeakage: 1
                    Lo2IF_dBc: 30
                    Lo2RF_dBc: 30
                         fNCO: 0
              flag_IQimb_comp: 0
               imb_AmpDB_comp: 0
              imb_PhsDeg_comp: 0
        flag_IQLvlOffset_comp: 0
```

**Waveform Upconversion**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/1b1a576c-c4c0-45b9-a85c-607ea8012522)

**Receiver Waveform Demodulation EVM Result**
- TM3.1
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/930d7343-b456-468d-a84a-4c42b4b9423f)
- TM1.1
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/5c8f3150-acd6-488b-93f1-6ad6c0f9049e)

**Receiver Waveform Demodulation Constellation Result**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/ae307ac8-eaab-4b2d-9587-815f6f74dc30)

**Receiver Waveform CCDF Result**
![image](https://github.com/kaycelin/ica_RF_Lo_Mixer_QEC/assets/87049112/3900bb03-9c3b-4ebf-9690-5c4e0d81de70)












