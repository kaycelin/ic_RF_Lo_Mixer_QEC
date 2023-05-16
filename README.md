# 0302_RF_Lo_Mixer_QEC

## Design flow
![image](https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/fd3ff425-95fb-4302-a440-0fee36630a63)


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
                     flag_IQimbalance: 1
                            imb_AmpDB: -3
                           imb_PhsDeg: 15
                     flag_IQLvlOffset: 1
                      imb_LvlOffsetDB: 3.1
                      flag_IQimb_comp: 0
                       imb_AmpDB_comp: 0
                      imb_PhsDeg_comp: 0
                flag_IQLvlOffset_comp: 0
                                 fnum: 41901
                                   fs: 1966080000
                               pwrdBm: 75.5631538720216
                                 fNCO: 100000000
                           IFfilterBW: []
                          IFfilterTyp: []
      ```
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/e44d7744-2ab7-43dc-90a6-ea2719cfcc47" width="60%">

  - **Downconversion**
    - Generate Lo for downconversion
      ```js
        LoClsD = 

          Lo with properties:

                            fs: 1966080000
                             T: 0.001
                           fLo: 491520000
                        Nsamps: 983040
                   flag_LoQuad: 'IQ'
            flag_LoQuad_method: 'LoIQ->LoPN'
             flag_LoPhaseNoise: 0
                        PwrdBm: 0
             phsNzFreqOffsetHz: [100000 200000 400000 600000 800000 1200000 … ]
                   phsNzLvldBc: [-90 -95 -100 -105 -110 -115 -120 -125 -130]
                             y: [983040×2 double]
                            yI: []
                            yQ: []
                          fnum: 41802
      ```
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/73b9f1fa-9c49-4c62-b250-1d16bdbe46b5" width="80%">
    - Down-mixer parameters and setting
      ```js
        if flag_imbTxEstimation
            % set digital mixer fNCO = 0 to estimate the iq imbalance of amp. and phase
            mixParmsD.fNCO = 0;
            
            % apply IF filter to supress mixer's RF image for the iq imbalance estimation
            mixParmsD.IFfilterBW = 0.5*(fs/2-100e6)*[-1 1]
        else
            mixParmsD.fNCO = -mixParms.fNCO;

            % apply IF filter to supress mixer's RF image for the demodulation
            mixParmsD.IFfilterBW = 0.5*(fs/2-100e6)*[-1 1];
        end
      ```
    - Downconversion: 
      ```js
        mixParmsD = 

          struct with fields:

                                    x: [983040×1 double]
                                   lo: [983040×2 double]
                                   fs: 1966080000
                                 fnum: 41902
                       flag_UDconvert: 'D'
            flag_UDconvert_FreqSelect: 'H'
                         flag_MixType: 'IQ'
                     flag_IQimbalance: 0
                            imb_AmpDB: -3
                           imb_PhsDeg: 15
                     flag_IQLvlOffset: 0
                      imb_LvlOffsetDB: 3.1
                             flag_Amp: 0
                            flag_Amp2: 0
                        mixerOut_flat: 0
                       flag_LoLeakage: 0
                            Lo2IF_dBc: 30
                            Lo2RF_dBc: 30
                                 fNCO: 0
                      flag_IQimb_comp: 0
                       imb_AmpDB_comp: 0
                      imb_PhsDeg_comp: 0
                flag_IQLvlOffset_comp: 0
                           IFfilterBW: [-441520000 441520000]
      ```
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/bfadd17b-77bb-4e9e-8c67-9b669aa0af52" width="60%">
  
  - **Tx Imbalance Calculation**
      ```js
        [IMB_MagDB_TXest, IMB_PhsDeg_TXest] = IQ_IMB_cor(real(yIfFir),  imag(yIfFir), 'TXQEC')
      ```
      - Result
      ```js
        IMB_MagDB_TXest = -2.88964490603545
        IMB_PhsDeg_TXest = 15.8632359822716
      ``` 
  - **DSP and Demodulation**
      ```js
        % NCO
        yNCO = nco(yADCin, -mixParms.fNCO, fs, [], 051002);
        
        % channel filter
        bwCh = (50e6/2+50e6)*[-1 1]
        [bChfir yNcoCh] = firGen(yNCO, fs, bwCh, 1e6, 'lpf', 'eqrip', [0.001 50], 051003);
        
        % evm check
        evmYNcoCh = evmInband(yNCO, x, fs, bwInband)
        evmYNcoChFir = evmInband(yNcoCh, x, fs, [])
      ```
      - Result
      ```js
        evmYNcoCh = 2.25509498307382
        evmYNcoChFir = 2.25677497715655
      ``` 
  - **Demodulation and EVM**
      ```js
        flagEVM = 'Matlab' % Input

        switch flagEVM % Demodulation and EVM
            case {'RAW','Raw','raw'}
                evmYDemod = evmInband(rxWaveform, x, fs, bwInband)

            case {'MATLAB','Matlab','matlab'}
                if 1 % Initial
                    evm3GPP = false; % |evm3GPP| is disabled for PDCCH EVM measurement.
                    targetRNTIs = []; % The example calculates the PDSCH EVM for the RNTIs listed above. To override the default PDSCH RNTIs, specify the |targetRNTIs| vector
                    plotEVM = true;
                    displayEVM = true;
                end
                cfg = struct();
                cfg.Evm3GPP = evm3GPP;
                cfg.TargetRNTIs = targetRNTIs;
                cfg.PlotEVM = plotEVM;
                cfg.DisplayEVM = displayEVM;
                cfg.Label = dlrefwavegen.ConfiguredModel{1};
                Ncarriers = Info.Ncarriers;
                for k = 1:Ncarriers
                    % Compute and display EVM measurements
                    [evmInfo,eqSym,refSym] = hNRDownlinkEVM_k(dlrefwavegen.Config,rxWaveform(:,k),cfg);
                    evm_PDCCH_RMS(k) = 100*evmInfo.PDCCH.OverallEVM.RMS
                    evm_PDSCH_RMS(k) = 100*evmInfo.PDSCH.OverallEVM.RMS
                end
            case {'kc'}
                dlRefwaveConfig = dlrefwavegen.Config;
                dlRefwaveInfo = dlrefwaveinfo;
                dlResourceInfo = dlresourceinfo;
                evmNRdl(dlRefwaveConfig, dlRefwaveInfo, dlResourceInfo, rxWaveform);
        end
      ```
      - Result
      ```js
        evm_PDCCH_RMS =  0.162929129123985
        evm_PDSCH_RMS = 1.34080953343171
      ``` 
      <img src="https://github.com/kaycelin/0302_RF_Lo_Mixer_QEC/assets/87049112/68430782-011b-4bc8-8828-fd6879889200" width="70%">

