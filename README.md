# 0302_RF_Lo_Mixer_QEC

## Design flow
<img src="https://user-images.githubusercontent.com/87049112/226227322-193bd0cb-974c-423a-a995-c5f428c9c972.png" width="80%">
- Summary_NR Test model.md


## Main: Mixer_Lo_main_k0510.m
  - **Lo and Phase Noise**
    - Lo parameters and settings
      ```js
      LoParms.phNzFreqOff: Phase noise freq. offset [Hz]
      LoParms.phNzLevl: Phase noise level [dBm/Hz]
      LoParms.fs
      LoParms.fLo
      LoParms.Nsamps
      LoParms.flag_LoQuad_method: 'LoPN->LoIQ' or 'LoIQ->LoPN'
      LoParms.PwrdBm: Sine wave power level [dBm]
      LoParms.flag_LoQuad: 'Real' (0 degree phase diff.) or 'IQ' (90 degree phase diff.)
      flag_LoPhaseNoiseModel: 'Default' or 'NR Phase Noise Model'
      ```
    - Generate Lo 
      ```js
        Lo with properties:

                        fs: 1966080000
                         T: 0.001
                       fLo: 491520000
                    Nsamps: 983040
               flag_LoQuad: 'IQ'
        flag_LoQuad_method: 'LoIQ->LoPN'
         flag_LoPhaseNoise: 1
                    PwrdBm: 0
         phsNzFreqOffsetHz: [100000 200000 400000 600000 800000 1200000 1800000 6000000 10000000]
               phsNzLvldBc: [-90 -95 -100 -105 -110 -115 -120 -125 -130]
                         y: [983040×2 double]
                        yI: []
                        yQ: []
                      fnum: 41802
      ```
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/4aaba4c0-ad5f-4ebb-968d-495a0f0bc816" width="80%">

  - **Upconversion**
    - Up-mixer parameters and setting
      ```js
        mixParms.x: Input signal
        mixParms.lo: Input Lo
        mixParms.fs
        ======= mixer type/iq imbalance/dc offset ========
        mixParms.flag_UDconvert: 'U' or 'D'
        mixParms.flag_UDconvert_FreqSelect: 'H' or 'L'
        mixParms.flag_MixType: 'IQ' or 'Real'/'IF'
        mixParms.flag_IQimbalance: 0 or 1
        mixParms.imb_AmpDB: Amplitude imbalance [dB]
        mixParms.imb_PhsDeg: Phase imbalance [degree]
        mixParms.flag_IQLvlOffset: 0 or 1, IQ level offset
        mixParms.imb_LvlOffsetDB: IQ mixer level offset [dB]
        ====== mixer unlinearity or amplifier stage ========
        mixParms.flag_Amp: Mixer unlinearity
        mixParms.ampIn_op1dB_dBm
        mixParms.ampIn_oip3_dBm
        mixParms.ampIn_gain_dB
        mixParms.ampOut_nf_dB
        mixParms.ampOut_flat: [freqList, flatnessList], [Hz, dB]
        mixParms.mixerOut_flat: [freqList, flatnessList], [Hz, dB]
        == mixer unlinearity or amplifier stage, mixer output ==
        mixParms.flag_Amp2: Mixer unlinearity
        mixParms.amp2In_op1dB_dBm
        mixParms.amp2In_oip3_dBm
        mixParms.amp2In_gain_dB
        mixParms.amp2Out_nf_dB
        mixParms.amp2Out_flat: [freqList, flatnessList], [Hz, dB]
        ============ mixer Lo leakage ===============
        mixParms.flag_LoLeakage: 0 or 1
        mixParms.Lo2IF_dBc: Lo leakage to IF [dBc]
        mixParms.Lo2RF_dBc: Lo leakage to RF [dBc]
        ============= others settings ================
        mixParms.fNCO: digital mixer [Hz]
        mixParms.flag_IQimb_comp: pre-compensate iq imbalance
        mixParms.imb_AmpDB_comp
        mixParms.imb_PhsDeg_comp
      ```
    - Upconversion
      ```js
        mixerCls = 

        mixer with properties:

                                x: [983040×1 double]
                               lo: [983040×2 double]
                                y: [983040×1 double]
                               yI: [983040×1 double]
                               yQ: [983040×1 double]
                   flag_UDconvert: 'U'
        flag_UDconvert_FreqSelect: 'H'
                     flag_MixType: 'IQ'
                         flag_Amp: 0
                  ampIn_op1dB_dBm: 1000
                   ampIn_oip3_dBm: 1000
                    ampIn_gain_dB: 0
                     ampOut_nf_dB: 0
                      ampOut_flat: []
                        flag_Amp2: 0
                 amp2In_op1dB_dBm: 1000
                  amp2In_oip3_dBm: 1000
                   amp2In_gain_dB: 0
                    amp2Out_nf_dB: 0
                     amp2Out_flat: []
                    mixerOut_flat: []
                   flag_LoLeakage: 1
                        Lo2IF_dBc: 30
                        Lo2RF_dBc: 30
                 flag_IQimbalance: 0
                        imb_AmpDB: -2
                       imb_PhsDeg: 0.5
                 flag_IQLvlOffset: 1
                  imb_LvlOffsetDB: 3.1
                  flag_IQimb_comp: 0
                   imb_AmpDB_comp: 0
                  imb_PhsDeg_comp: 0
            flag_IQLvlOffset_comp: 0
                             fnum: 41901
                               fs: 1966080000
                           pwrdBm: 75.3091071205401
                             fNCO: 0
                       IFfilterBW: []
                      IFfilterTyp: []
      ```
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/b2c12d51-8e12-468a-98a1-f7491f178260" width="80%">

  - **Downconversion**
    - Down-mixer parameters and setting
      ```js
        mixParmsD.fNCO: 0
        mixParms.lo: Input Lo
        mixParms.fs
        ======= mixer type/iq imbalance/dc offset ========
        mixParms.flag_UDconvert: 'U' or 'D'
        mixParms.flag_UDconvert_FreqSelect: 'H' or 'L'
        mixParms.flag_MixType: 'IQ' or 'Real'/'IF'
        mixParms.flag_IQimbalance: 0 or 1
        mixParms.imb_AmpDB: Amplitude imbalance [dB]
        mixParms.imb_PhsDeg: Phase imbalance [degree]
        mixParms.flag_IQLvlOffset: 0 or 1, IQ level offset
        mixParms.imb_LvlOffsetDB: IQ mixer level offset [dB]
        ====== mixer unlinearity or amplifier stage ========
        mixParms.flag_Amp: Mixer unlinearity
        mixParms.ampIn_op1dB_dBm
        mixParms.ampIn_oip3_dBm
        mixParms.ampIn_gain_dB
        mixParms.ampOut_nf_dB
        mixParms.ampOut_flat: [freqList, flatnessList], [Hz, dB]
        mixParms.mixerOut_flat: [freqList, flatnessList], [Hz, dB]
        == mixer unlinearity or amplifier stage, mixer output ==
        mixParms.flag_Amp2: Mixer unlinearity
        mixParms.amp2In_op1dB_dBm
        mixParms.amp2In_oip3_dBm
        mixParms.amp2In_gain_dB
        mixParms.amp2Out_nf_dB
        mixParms.amp2Out_flat: [freqList, flatnessList], [Hz, dB]
        ============ mixer Lo leakage ===============
        mixParms.flag_LoLeakage: 0 or 1
        mixParms.Lo2IF_dBc: Lo leakage to IF [dBc]
        mixParms.Lo2RF_dBc: Lo leakage to RF [dBc]
        ============= others settings ================
        mixParms.fNCO: digital mixer [Hz]
        mixParms.flag_IQimb_comp: pre-compensate iq imbalance
        mixParms.imb_AmpDB_comp
        mixParms.imb_PhsDeg_comp
      ```
