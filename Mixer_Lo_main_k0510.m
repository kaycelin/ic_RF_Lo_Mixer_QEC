% *Input Settings*
%% 
% * flag_SigType: 'nr' or 'cw'
% * flag_LoPhaseNoise: 0 or 1
% * =============================================
% * fs: Sampling freq.
% * fLo: Local oscillator freq.

flag_SigType = 'nr'
flag_LoPhaseNoise = 1

% *Signal*

% Input signal
try
    switch flag_SigType
        case 'nr'
            fs = fs;
            Nsamps = numel(x);
    end
catch
    load('100MHz_scs30_fs983p04MHz_TM31_cfr8p5_1frame.mat')
    fs = Info.fs;
    bwInband = Info.bwInband;
    dlresourceinfo = Info;
    dlrefwavegen = Info.dlrefwavegen;
    dlresourceinfo = Info.dlresourceinfo;
    x = IQ;
    Nsamps = numel(x);

end
PLOT_FFT_dB(x, fs, [], 'dBm', 041801, {[],flag_SigType})

fLo = fs/4

% *Lo and Phase Noise*
%% 
% * LoParms.phNzFreqOff: Phase noise freq. offset [Hz]
% * LoParms.phNzLevl: Phase noise level [dBm/Hz]
% * LoParms.fs
% * LoParms.fLo
% * LoParms.Nsamps
% * LoParms.flag_LoQuad_method: 'LoPN->LoIQ' or 'LoIQ->LoPN'
% * LoParms.PwrdBm: Sine wave power level [dBm]
% * LoParms.flag_LoQuad: 'Real' (0 degree phase diff.) or 'IQ' (90 degree phase 
% diff.)
% * ==============================================================
% * flag_LoPhaseNoiseModel: 'Default' or 'NR Phase Noise Model'

% Input of Lo parameters
LoParms.flag_LoPhaseNoise = flag_LoPhaseNoise;
if 1
    flag_LoPhaseNoiseModel = 'Default'
else
    flag_LoPhaseNoiseModel = 'NR Phase Noise Model'
end

if LoParms.flag_LoPhaseNoise
    switch flag_LoPhaseNoiseModel
        case 'Default'
            LoParms.flag_LoPhaseNoise = 1;

            deltaPN = 1
            if LoParms.flag_LoPhaseNoise
                phNzFreqOff = [ 100e3, 200e3, 400e3, 600e3, 800e3,...
                    1.2e6, 1.8e6, 6e6, 10e6 ]; % Offset From Carrier
                phNzLevel = [ -100, -105, -110, -115, -120,...
                    -125, -130, -135, -140 ] + deltaPN*10; % Phase Noise power

                LoParms.phsNzFreqOffsetHz = phNzFreqOff;
                LoParms.phsNzLvldBc = phNzLevel;
            end
        case 'NR Phase Noise Model'
            LoParms.flag_LoPhaseNoise = 0;

            phNzFreqOff = (4:0.2:log10(fs/2)); % Model offset from 1e4 Hz to sr/2 Hz
            foffset = 10.^phNzFreqOff;         % Linear frequency offset
            PNModel = 'A'; % 'A' (TDoc R1-163984 Set A), 'B' (TDoc R1-163984 Set B), 'C' (TR 38.803)
            pn_PSD = hPhaseNoisePoleZeroModel(foffset,fLo,PNModel); % dBc/Hz
            % Set phase noise level
            pnoise = comm.PhaseNoise('FrequencyOffset',foffset,'Level',pn_PSD,'SampleRate',fs);
            pnoise.RandomStream = "mt19937ar with seed";

            if 1
                % Visualize spectrum mask of phase noise
                figure(041803,1,2,1)
                semilogx(foffset,pn_PSD)
                xlabel('Frequency offset (Hz)')
                ylabel('dBc/Hz')
                title('Phase noise magnitude response')
                grid on
            end
    end
end

LoParms.fs = fs;
LoParms.fLo = fLo;
LoParms.Nsamps = Nsamps;
LoParms.flag_LoQuad_method = 'LoIQ->LoPN';
LoParms.PwrdBm = 0;
LoParms.fnum = 041802;
if 0
    LoParms.flag_LoQuad = 'Real';
else
    LoParms.flag_LoQuad = 'IQ';
end

% Generate Lo
[LoCls,LoStruct] = Lo(LoParms);

% Export Lo
LoSignal = LoCls.y;

% Apply Phase Noise Model
switch flag_LoPhaseNoiseModel
    case 'NR Phase Noise Model'
        LoSignal = pnoise(LoSignal);

        % Plot
        PLOT_FFT_dB(LoSignal(:,1),fs,[],'dBc',[041803,1,2,2],'LoSignal',[],...
            'RBW',[],1,1e6,0,[fLo, foffset+fLo]);
end

% *UpDownConvert Settings*
%% 
% * mixParms.x: Input signal
% * mixParms.lo: Input Lo
% * mixParms.fs
% * ======= mixer type/iq imbalance/dc offset ========
% * mixParms.flag_UDconvert: 'U' or 'D'
% * mixParms.flag_UDconvert_FreqSelect: 'H' or 'L'
% * mixParms.flag_MixType: 'IQ' or 'Real'/'IF'
% * mixParms.flag_IQimbalance: 0 or 1
% * mixParms.imb_AmpDB: Amplitude imbalance [dB]
% * mixParms.imb_PhsDeg: Phase imbalance [degree]
% * mixParms.flag_IQLvlOffset: 0 or 1, IQ level offset
% * mixParms.imb_LvlOffsetDB: IQ mixer level offset [dB]
% * ====== mixer unlinearity or amplifier stage ========
% * mixParms.flag_Amp: Mixer unlinearity
% * mixParms.ampIn_op1dB_dBm
% * mixParms.ampIn_oip3_dBm
% * mixParms.ampIn_gain_dB
% * mixParms.ampOut_nf_dB
% * mixParms.ampOut_flat: [freqList, flatnessList], [Hz, dB]
% * mixParms.mixerOut_flat: [freqList, flatnessList], [Hz, dB]
% * == mixer unlinearity or amplifier stage, mixer output ==
% * mixParms.flag_Amp2: Mixer unlinearity
% * mixParms.amp2In_op1dB_dBm
% * mixParms.amp2In_oip3_dBm
% * mixParms.amp2In_gain_dB
% * mixParms.amp2Out_nf_dB
% * mixParms.amp2Out_flat: [freqList, flatnessList], [Hz, dB]
% * ============ mixer Lo leakage ===============
% * mixParms.flag_LoLeakage: 0 or 1
% * mixParms.Lo2IF_dBc: Lo leakage to IF [dBc]
% * mixParms.Lo2RF_dBc: Lo leakage to RF [dBc]
% * ============= others settings ================
% * mixParms.fNCO: digital mixer [Hz]
% * mixParms.flag_IQimb_comp: pre-compensate iq imbalance
% * mixParms.imb_AmpDB_comp
% * mixParms.imb_PhsDeg_comp
% *Up Conversion*

% mixer signal and Lo
mixParms.x = x(:);
mixParms.lo = LoSignal;
mixParms.fs = fs;
mixParms.fnum = 041901

% mixer type/iq imbalance/dc offset
mixParms.flag_UDconvert = 'U';
mixParms.flag_UDconvert_FreqSelect = 'H';
if 0
    mixParms.flag_MixType = 'IQ';
else
    mixParms.flag_MixType = LoParms.flag_LoQuad;
end

mixParms.flag_IQimbalance = 1;
if mixParms.flag_IQimbalance
    mixParms.imb_AmpDB = -3*1;
    mixParms.imb_PhsDeg = 15*1;
end

mixParms.flag_IQLvlOffset = 1;
if mixParms.flag_IQLvlOffset
    mixParms.imb_LvlOffsetDB = 1*3.1;
end

% mixer unlinearity or amplifier stage
mixParms.flag_Amp = 0;
if mixParms.flag_Amp
    mixParms.ampIn_op1dB_dBm = 35;
    mixParms.ampIn_oip3_dBm = 5;
    mixParms.ampIn_gain_dB = 10;
    mixParms.ampOut_nf_dB = 5;
    mixParms.ampOut_flat = 0;
end

% mixer unlinearity or amplifier stage, mixer output
mixParms.flag_Amp2 = 0;
if mixParms.flag_Amp2
    mixParms.amp2In_op1dB_dBm = 1000;
    mixParms.amp2In_oip3_dBm = 1000;
    mixParms.amp2In_gain_dB = 0;
    mixParms.amp2Out_nf_dB = 0;
    if 1 % mixer flatness
        flatFreq = [0, fs/4, fs/2]
        flatLvlDB = [-10 0 -10]
        mixParms.amp2Out_flat = [flatFreq(:), flatLvlDB(:)];
    else
        mixParms.amp2Out_flat = [];
    end
end

if 0 % mixer flatness
    flatFreq = [0, fs/4, fs/2]
    flatLvlDB = [-2 0 -2]
    mixParms.mixerOut_flat = [flatFreq(:), flatLvlDB(:)];
else
    mixParms.mixerOut_flat = [];
end

% mixer Lo leakage
mixParms.flag_LoLeakage = 1;
if mixParms.flag_LoLeakage
    mixParms.Lo2IF_dBc = 30;
    mixParms.Lo2RF_dBc = 30;
end

% digital mixer, NCO
if 1
    mixParms.fNCO = 100e6;
else
    mixParms.fNCO = 0;
end

if 1 % pre-compensate iq imbalance
    if exist('IMB_MagDB_TXest', 'var')&&~isempty(IMB_MagDB_TXest)||...
            exist('IMB_PhsDeg_TXest', 'var')&&~isempty(IMB_PhsDeg_TXest)
        mixParms.flag_IQimb_comp = 1;
    else
        mixParms.flag_IQimb_comp = 0;
    end
    if mixParms.flag_IQimb_comp*mixParms.flag_IQimbalance
        if 0
            mixParms.imb_AmpDB_comp = + mixParms.imb_AmpDB;
            mixParms.imb_PhsDeg_comp = + mixParms.imb_PhsDeg;
        elseif exist('IMB_MagDB_TXest', 'var')&&~isempty(IMB_MagDB_TXest)||...
                exist('IMB_PhsDeg_TXest', 'var')&&~isempty(IMB_PhsDeg_TXest)
            mixParms.imb_AmpDB_comp = + IMB_MagDB_TXest;
            mixParms.imb_PhsDeg_comp = + IMB_PhsDeg_TXest;
        else
            mixParms.imb_AmpDB_comp = 0;
            mixParms.imb_PhsDeg_comp = 0;
        end
    else
        mixParms.imb_AmpDB_comp = 0;
        mixParms.imb_PhsDeg_comp = 0;
    end
end

if 1 % compensate iq level offset -> Not support here !
    mixParms.flag_IQLvlOffset_comp = 0;
end

% Generate Mixer class and mixing
[mixerCls,mixerStruct] = mixer(mixParms);

% Export
yRF = mixerCls.y;
PLOT_FFT_dB(yRF,fs,[],'dBm',1419,{'yRF','Upconversion'},[],[],[],[],1e6,1,[]);

% *Down Conversion, IF architecture, From RF down to IF*

% Initialization Mixer parameters
mixParmsD = mixParms;
if 1
    mixParmsD.x = yRF;
    if 0
        fIF = 100e6;
    else
        fIF = 0;
    end
    mixParmsD.flag_UDconvert = 'D';
    mixParmsD.flag_MixType = mixParms.flag_MixType;
    mixParmsD.flag_UDconvert_FreqSelect = mixParms.flag_UDconvert_FreqSelect;
    mixParmsD.flag_IQimbalance = 0;
    mixParmsD.flag_IQimb_comp = 0;
    mixParmsD.flag_IQLvlOffset = 0;
    mixParmsD.flag_LoLeakage = 0;
    mixParmsD.fnum = 041902;
    mixParmsD.mixerOut_flat = 0;
end

% Generate LoD for downconversion
LoParmsD = LoParms;
LoParmsD.flag_LoPhaseNoise = 0;
if 1
    LoParmsD.flag_LoQuad = LoParms.flag_LoQuad;
else
    LoParmsD.flag_LoQuad = 'real';
end
LoParmsD.fLo = LoParms.fLo - fIF;
[LoClsD,LoStructD] = Lo(LoParmsD);
if 1
    mixParmsD.lo = LoClsD.y;
end
if 1 
    flag_imbTxEstimation = 1

    % set digital mixer fNCO = 0 to estimate the iq imbalance of amp. and phase
    mixParmsD.fNCO = 0;
    
    % apply IF filter to supress mixer's RF image for the iq imbalance estimation
    mixParmsD.IFfilterBW = 0.5*(fs/2-100e6)*[-1 1]
else
    flag_imbTxEstimation = 0

    % apply IF filter to supress mixer's RF image for the demodulation
    mixParmsD.IFfilterBW = 0.5*(fs/2-100e6)*[-1 1];
    mixParmsD.fNCO = -mixParms.fNCO;
end

% Mixing process, from RF to IF
[mixerClsD,mixerStructD] = mixer(mixParmsD);
yIF = mixerClsD.y;
PLOT_FFT_dB(yIF,fs,[],'dBm',1419,{'yIF','Downconversion'},[],[],[],[],1e6,1,[]);

% *Imbalance calculation and ADC*

% Tx IQ imbalance estimation
if flag_imbTxEstimation
    if mixerClsD.fNCO ~= 0
        error('Calculate Imbalance before NCO !')
    elseif isempty(mixerClsD.IFfilterBW)
        % apply IF filter to supress mixer's RF image for the iq imbalance
        % estimation, make sure inband ripple is small !
        bwIF = 0.5*(fs/2-100e6)*[-1 1]
        [bIfFir, yIfFir] = firGen(yIF, fs, bwIF, 1e6, 'lpf', 'eqrip', [0.001 50], 051001);
    else
        yIfFir = yIF;
    end
    [IMB_MagDB_TXest, IMB_PhsDeg_TXest] = IQ_IMB_cor(real(yIfFir),  imag(yIfFir), 'TXQEC')
end

if 1 % ADC digital signal process
    if 1
        yADCin = yIfFir;
    else
        yADCin = yIF;
    end

    % NCO
    yNCO = nco(yADCin, -mixParms.fNCO, fs, [], 051002);

    % channel filter
    bwCh = (50e6/2+50e6)*[-1 1]
    [bChfir yNcoCh] = firGen(yNCO, fs, bwCh, 1e6, 'lpf', 'eqrip', [0.001 50], 051003);

    % evm check
    evmYNcoCh = evmInband(yNCO, x, fs, bwInband)
    evmYNcoChFir = evmInband(yNcoCh, x, fs, [])

    % export
    rxWaveform = yNcoCh;
end
% *Demodulation and EVM*

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