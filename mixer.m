%% 2022-11-18, draft
%% 2023-02-07, change flag_Mixtype to flag_MixType
%% 2023-02-07, flag_LoLeakage = [] default, LoLeakage_dBm
%% 2023-03-22, Amplifier effect: unlinearity, gain and noise figure, 2nd stage
%% 2023-05-05, add fNCO
%% 2023-05-05, export yI and yQ for RxQEC(QEC receiver)
%% 2023-05-10, flag_IQimb_Comp, pre-compensate the imbalance of amp. and phase
%% 2023-05-10, add IF filter before ADC process(Decimation/NCO/...)
%% 2023-05-11, flag_IQLvlOffset_comp, pre-compensate the iq level offset

classdef mixer < handle
    properties
        x = [] % input waveform
        lo = [] % lo waveform
        y = [] % output waveform
        yI = [] % 2023-05-05
        yQ = [] % 2023-05-05

        flag_UDconvert = 'U' 
        % U: upconversion
        % D: downconversion

        flag_UDconvert_FreqSelect = 'H'
        % H: high side freq.
        % L: low side freq.
        % Zero: zero freq. for downconversion

        flag_MixType = 'IQ'
        % IQ: iq mixer
        % Real: real mixer

        flag_Amp = 0 % Amplifier effect: unlinearity, gain and noise figure
        ampIn_op1dB_dBm = 1000
        ampIn_oip3_dBm = 1000
        ampIn_gain_dB = 0
        ampOut_nf_dB = 0
        ampOut_flat = []

        flag_Amp2 = 0 % Amplifier effect: unlinearity, gain and noise figure, 2nd stage, 2023-03-22
        amp2In_op1dB_dBm = 1000
        amp2In_oip3_dBm = 1000
        amp2In_gain_dB = 0
        amp2Out_nf_dB = 0
        amp2Out_flat = []

        mixerOut_flat = []

        flag_LoLeakage = 0
        Lo2IF_dBc = 0
        Lo2RF_dBc = 0

        flag_IQimbalance = 0
        imb_AmpDB = -2
        imb_PhsDeg = 0.5

        flag_IQLvlOffset = 0
        imb_LvlOffsetDB = 0

        flag_IQimb_comp = 0 %% 2023-05-10, pre-compensate the imbalance of amp. and phase
        imb_AmpDB_comp = 0
        imb_PhsDeg_comp = 0

        flag_IQLvlOffset_comp = 0 % 2023-05-11, pre-compensate the iq level offset

        fnum = [] % plot
        fs = [] % sampling rate

        pwrdBm = [] % power

        fNCO = 0 % digital mixer, 2023-05-05

        IFfilterBW = [] % 2023-05-10, add IF filter before ADC process(Decimation/NCO/...)
        IFfilterTyp = []

    end

    properties (Dependent)
    end

    methods (Access = public)
        function [obj, objStruct] = mixer(varargin)
            % Initialization and set parameters
            if iscell(varargin) && size(varargin,2)>1
                varargin_tmp = varargin;
            else
                varargin_tmp = varargin{:};
            end
            obj.setParameters(varargin_tmp);

            % mixer unlinearity
            if obj.flag_Amp
                lna = [];
                lna.model = "Cubic polynomial"; % model of am/am and am/pm
                lna.gain_dB = obj.ampIn_gain_dB; % Input gain_dB
                lna.oip3_dBm = obj.ampIn_oip3_dBm; % Input oip3
                lna.op1dB_dBm = obj.ampIn_op1dB_dBm;
                lna.opsat_dBm = lna.op1dB_dBm;
                lna.TOISpecification = 'OP1dB';
                lna.lna_NFdB = obj.ampOut_nf_dB;
                lna.flat_dB = obj.ampOut_flat;
                lna.flag_lna = 'LNA';
                if isempty(obj.y)
                    lna.x = obj.x;
                else
                    lna.x = obj.y;
                end
                lna.fs = obj.fs;
                lna.plt = [];
                if 1
                    lna.plt.fnum = [];
                end
                [lnaCls,~] = powerAmp(lna);
                obj.y = lnaCls.y;
            end

            % mixing
            obj.mix();
            
            % mixer unlinearity or flatness
            if obj.flag_Amp2 || ~isempty(obj.mixerOut_flat)*(~isscalar(obj.mixerOut_flat))
                lna2 = [];
                lna2.model = "Cubic polynomial"; % model of am/am and am/pm
                lna2.gain_dB = obj.amp2In_gain_dB; % Input gain_dB
                lna2.oip3_dBm = obj.amp2In_oip3_dBm; % Input oip3
                lna2.op1dB_dBm = obj.amp2In_op1dB_dBm;
                lna2.opsat_dBm = lna2.op1dB_dBm;
                lna2.TOISpecification = 'OP1dB';
                lna2.lna_NFdB = obj.amp2Out_nf_dB;
                try
                    lna2.flat_dB = obj.mixerOut_flat;
                catch
                    lna2.flat_dB = obj.amp2Out_flat;
                end
                lna2.flag_lna = 'LNA';
                if isempty(obj.y)
                    lna2.x = obj.x;
                else
                    lna2.x = obj.y;
                end
                lna2.fs = obj.fs;
                lna2.plt = [];
                if 1
                    lna2.plt.fnum = [];
                end
                [lna2Cls,~] = powerAmp(lna2);
                obj.y = lna2Cls.y;
            end

            if 0 % flatness, flatness done on lna2 stage
                if ~isempty(obj.amp_flat)
                    flatdB_In = obj.amp_flat;
                    [y_flat] = lnaCls.pa_flat(obj.y, flatdB_In, obj.fs);
                    obj.y = y_flat;
                    % plot
                    if ~isempty(obj.fnum)
                        obj.plot(obj.y, obj.fs, 'dBm', 'kHz', 0, 'flatness', obj.fnum);
                    end
                end
            end

            % plot
            if ~isempty(obj.fnum)
                obj.plot(obj.y, obj.fs, 'dBm', 'kHz', 0, {'mixer output','Mixer'}, obj.fnum);
            end
            % export to structure
            objStruct = struct(obj);

        end
    end

    properties (Access = private)
        % flag_setParameters = 0, reset for setting
        flag_setParameters = 0
        default
        Nstages = 1
        typeInputSignal = 'real' % 'complx' 
    end

    methods
        function setParameters(obj, Inputs, n1)
            warning off MATLAB:structOnObject
            obj.default = struct(obj);

            % format type
            if isstruct(Inputs)
                InputTmp = Inputs;
                nEnd = 1;
            elseif iscell(Inputs) && size(Inputs, 1)>1
                InputTmp = Inputs;
                nEnd = size(InputTmp,1);
            elseif iscell(Inputs)
                InputTmp = Inputs;
                nEnd = size(InputTmp,1);
            else
                error('mixer.setParameters: type of params is STRUCT or CELL!')
            end

            % multi stages
            if ~exist('n1', 'var')||isempty(n1)
                n1 = 1;
                obj.Nstages = nEnd;
                flag_setParams = 1;
            else
                nEnd = n1;
                flag_setParams = 0;
            end

            % multi stages loop
            for n=n1:nEnd
                % multi conditions parameters
                Nvar = size(InputTmp,2);

                % format type: params_tmp, format transfer: varargin_tmp
                i = 1;
                if iscell(InputTmp)
                    mixParams = InputTmp{n,i};
                elseif isstruct(InputTmp)
                    mixParams = InputTmp;
                else
                    error('powerAmp, setParameters: type of params type is STRUCT or CELL!')
                end

                if Nvar>i-1 && ~isempty(mixParams)
                    if isstruct(mixParams)
                        % mixer parameters
                        if isfield(mixParams,'x')&&~isempty(mixParams.x)
                            if size(mixParams.x,2)>size(mixParams.x,1)
                                obj.x = mixParams.x.';
                            else
                                obj.x = mixParams.x;
                            end
                        else
                            error('mixer.setParameters: Input signal x should NOT be EMPTY!')
                        end
                        if isfield(mixParams,'lo')&&~isempty(mixParams.lo)
                            obj.lo = mixParams.lo;
                        end
                        if 0
                            if isfield(mixParams,'y')&&~isempty(mixParams.y)
                                obj.y = mixParams.y;
                            end
                        end
                        if isfield(mixParams,'flag_UDconvert')&&~isempty(mixParams.flag_UDconvert)
                            obj.flag_UDconvert = mixParams.flag_UDconvert;
                            switch obj.flag_UDconvert
                                case 'U'
                                    obj.flag_UDconvert_FreqSelect = 'H';
                                case 'D'
                                    obj.flag_UDconvert_FreqSelect = 'L';
                            end
                        end
                        if isfield(mixParams,'flag_UDconvert_FreqSelect')&&~isempty(mixParams.flag_UDconvert_FreqSelect)
                            obj.flag_UDconvert_FreqSelect = mixParams.flag_UDconvert_FreqSelect;
                        end
                        if isfield(mixParams,'flag_MixType')&&~isempty(mixParams.flag_MixType)
                            obj.flag_MixType = mixParams.flag_MixType;
                        end
                        if isfield(mixParams,'flag_IQimbalance')&&~isempty(mixParams.flag_IQimbalance)
                            obj.flag_IQimbalance = mixParams.flag_IQimbalance;
                        end
                        if isfield(mixParams,'imb_AmpDB')&&~isempty(mixParams.imb_AmpDB)
                            obj.imb_AmpDB = mixParams.imb_AmpDB;
                        end
                        if isfield(mixParams,'imb_PhsDeg')&&~isempty(mixParams.imb_PhsDeg)
                            obj.imb_PhsDeg = mixParams.imb_PhsDeg;
                        end
                        if isfield(mixParams,'flag_IQLvlOffset')&&~isempty(mixParams.flag_IQLvlOffset)
                            obj.flag_IQLvlOffset = mixParams.flag_IQLvlOffset;
                        end
                        if isfield(mixParams,'imb_LvlOffsetDB')&&~isempty(mixParams.imb_LvlOffsetDB)
                            obj.imb_LvlOffsetDB = mixParams.imb_LvlOffsetDB;
                        end
                        if isfield(mixParams,'fnum')&&~isempty(mixParams.fnum)
                            obj.fnum = mixParams.fnum;
                        end
                        if isfield(mixParams,'fs')&&~isempty(mixParams.fs)
                            obj.fs = mixParams.fs;
                        end
                        % 2023-02-07
                        if isfield(mixParams,'flag_LoLeakage')&&~isempty(mixParams.flag_LoLeakage)
                            obj.flag_LoLeakage = mixParams.flag_LoLeakage;
                        end
                        if isfield(mixParams,'Lo2IF_dBc')&&~isempty(mixParams.Lo2IF_dBc)
                            obj.Lo2IF_dBc = mixParams.Lo2IF_dBc;
                        end
                        if isfield(mixParams,'Lo2RF_dBc')&&~isempty(mixParams.Lo2RF_dBc)
                            obj.Lo2RF_dBc = mixParams.Lo2RF_dBc;
                        end
                        % 2023-02-20
                        if 1
                            if isfield(mixParams,'flag_Amp')&&~isempty(mixParams.flag_Amp)
                                obj.flag_Amp = mixParams.flag_Amp;
                            end
                            if isfield(mixParams,'ampIn_op1dB_dBm')&&~isempty(mixParams.ampIn_op1dB_dBm)
                                obj.ampIn_op1dB_dBm = mixParams.ampIn_op1dB_dBm;
                            end
                            if isfield(mixParams,'ampIn_oip3_dBm')&&~isempty(mixParams.ampIn_oip3_dBm)
                                obj.ampIn_oip3_dBm = mixParams.ampIn_oip3_dBm;
                            end
                            if isfield(mixParams,'ampIn_gain_dB')&&~isempty(mixParams.ampIn_gain_dB)
                                obj.ampIn_gain_dB = mixParams.ampIn_gain_dB;
                            end
                            if isfield(mixParams,'ampOut_nf_dB')&&~isempty(mixParams.ampOut_nf_dB)
                                obj.ampOut_nf_dB = mixParams.ampOut_nf_dB;
                            end
                            if isfield(mixParams,'ampOut_flat')&&~isempty(mixParams.ampOut_flat)
                                obj.ampOut_flat = mixParams.ampOut_flat;
                            end
                        end

                        if 1
                            if isfield(mixParams,'flag_Amp2')&&~isempty(mixParams.flag_Amp2)
                                obj.flag_Amp2 = mixParams.flag_Amp2;
                            end
                            if isfield(mixParams,'amp2In_op1dB_dBm')&&~isempty(mixParams.amp2In_op1dB_dBm)
                                obj.amp2In_op1dB_dBm = mixParams.amp2In_op1dB_dBm;
                            end
                            if isfield(mixParams,'amp2In_oip3_dBm')&&~isempty(mixParams.amp2In_oip3_dBm)
                                obj.amp2In_oip3_dBm = mixParams.amp2In_oip3_dBm;
                            end
                            if isfield(mixParams,'amp2In_gain_dB')&&~isempty(mixParams.amp2In_gain_dB)
                                obj.amp2In_gain_dB = mixParams.amp2In_gain_dB;
                            end
                            if isfield(mixParams,'amp2Out_nf_dB')&&~isempty(mixParams.amp2Out_nf_dB)
                                obj.ampOut_nf_dB = mixParams.ampOut_nf_dB;
                            end
                            if isfield(mixParams,'amp2Out_flat')&&~isempty(mixParams.amp2Out_flat)
                                obj.amp2Out_flat = mixParams.amp2Out_flat;
                            end
                            % mixer flatness
                            if isfield(mixParams,'mixerOut_flat')&&~isempty(mixParams.mixerOut_flat)
                                obj.mixerOut_flat = mixParams.mixerOut_flat;
                            end
                        end

                        if 1 % 2023-05-05
                            if isfield(mixParams,'fNCO')&&~isempty(mixParams.fNCO)
                                obj.fNCO = mixParams.fNCO;
                            end
                        end

                        if 1 % 2023-05-10, flag_IQimb_Comp, pre-compensate the imbalance of amp. and phase
                            if isfield(mixParams,'flag_IQimb_comp')&&~isempty(mixParams.flag_IQimb_comp)
                                obj.flag_IQimb_comp = mixParams.flag_IQimb_comp;
                            end
                            if isfield(mixParams,'imb_AmpDB_comp')&&~isempty(mixParams.imb_AmpDB_comp)
                                obj.imb_AmpDB_comp = mixParams.imb_AmpDB_comp;
                            end
                            if isfield(mixParams,'imb_PhsDeg_comp')&&~isempty(mixParams.imb_PhsDeg_comp)
                                obj.imb_PhsDeg_comp = mixParams.imb_PhsDeg_comp;
                            end

                            if isfield(mixParams,'flag_IQLvlOffset_comp')&&~isempty(mixParams.flag_IQLvlOffset_comp)
                                obj.flag_IQLvlOffset_comp = mixParams.flag_IQLvlOffset_comp;
                            end
                            
                        end

                        if 1 % 2023-05-10, add IF filter before ADC process(Decimation/NCO/...)
                            if isfield(mixParams,'IFfilterBW')&&~isempty(mixParams.IFfilterBW)
                                obj.IFfilterBW = mixParams.IFfilterBW;
                            end
                            if isfield(mixParams,'IFfilterTyp')&&~isempty(mixParams.IFfilterTyp)
                                obj.IFfilterTyp = mixParams.IFfilterTyp;
                            end        
                        end
                    else
                        error('mixer.setParameters: type of mixParams type should be STRUCT!')
                    end
                end
                i = i+1;
            end
        end
    end

    methods
        function [y] = mix(obj, x, lo, flag_UDconvert, flag_MixType, flag_IQimbalance, flag_LoLeakage, fnum)
            if ~exist('x', 'var') || isempty(x)
                x = obj.x;
            end

            if ~exist('lo', 'var') || isempty(lo)
                lo = obj.lo;
            end

            if ~exist('flag_UDconvert', 'var') || isempty(flag_UDconvert)
                flag_UDconvert = obj.flag_UDconvert;
            end

            if ~exist('flag_MixType', 'var') || isempty(flag_MixType)
                flag_MixType = obj.flag_MixType;
            end

            if strcmpi(flag_MixType, 'Real')
                flag_IQimbalance = 0;
                flag_IQLvlOffset = 0;
            end

            if ~exist('flag_IQimbalance', 'var') || isempty(flag_IQimbalance)
                flag_IQimbalance = obj.flag_IQimbalance;
            elseif numel(flag_IQimbalance)==2
                flag_IQimbalance = 1;
                obj.imb_AmpDB = flag_IQimbalance(1);
                obj.imb_PhsDeg = flag_IQimbalance(2);
            else
                flag_IQimbalance = 0;
            end

            if ~exist('flag_IQLvlOffset', 'var') || isempty(flag_IQLvlOffset)
                flag_IQLvlOffset = obj.flag_IQLvlOffset;
            elseif obj.imb_LvlOffsetDB == 0
                obj.imb_LvlOffsetDB = flag_IQLvlOffset;
            else
                flag_IQLvlOffset = 0;
            end

            if ~exist('flag_LoLeakage', 'var') || isempty(flag_LoLeakage)
                flag_LoLeakage = obj.flag_LoLeakage;
            elseif numel(flag_LoLeakage)==2
                flag_LoLeakage = 1;
                obj.Lo2IFdBc = flag_LoLeakage(1);
                obj.Lo2RFdBc = flag_LoLeakage(2);
            else
                flag_LoLeakage = 0;
            end

            if ~exist('fs', 'var') || isempty(fs)
                fs = obj.fs;
            end
            if ~exist('fnum', 'var') || isempty(obj.fs)
                fnum = [];
            end

            if ~isreal(x)
                if abs(pwrdB(imag(x))-pwrdB(real(x))) <= 20 % double check power of i and q
                    obj.typeInputSignal = 'complx';
                end
            end

            if obj.fNCO ~= 0 && strcmpi(obj.flag_UDconvert, 'U')
                x = nco(x, obj.fNCO, fs, [], []);
            end

            % column
            switch obj.flag_UDconvert
                case 'U'
                    if any(size([x])==1) && isrow(x)
                        flag_row = 1;
                        xIF = x(:);
                        xI = real(xIF);
                        xQ = imag(xIF);
                    elseif any(size([x])==1) && ~isrow(x)
                        flag_row = 0;
                        xIF = x;
                        xI = real(xIF);
                        xQ = imag(xIF);
                    elseif size(x,2)>size(x,1)
                        flag_row = 1;
                        xI = x(1,:).';
                        xQ = x(2,:).';
                        xIF = xI + 1i*xQ;
                    elseif size(x,2)<size(x,1)
                        flag_row = 0;
                        xI = x(:,1);
                        xQ = x(:,2);
                        xIF = xI + 1i*xQ;
                    end
                case 'D'
                    if size(x,2)>size(x,1)
                        flag_row = 1;
                    else
                        flag_row = 0;
                    end
                    xIF = x/sqrt(2);

            end

            % column, lo
            if any(size([lo])==1) 
                if 0
                    if ~strcmpi(flag_MixType, 'Real')
                        error('mixer.mix: check the lo input parameter ?!')
                    end
                end
                loI = lo;
                loQ = lo;
            elseif size(lo,2)>size(lo,1)
                loI = lo(1,:).';
                loQ = lo(2,:).';
            elseif size(lo,2)<size(lo,1)
                loI = lo(:,1);
                loQ = lo(:,2);
            end

            % lo leakage to IF
            if flag_LoLeakage
                Lo2IFdBc = abs(obj.Lo2IF_dBc);
                if Lo2IFdBc~=0
                    Lo2IF = lo / 10^(Lo2IFdBc/20);
                    Lo2IF_i = loI / 10^(Lo2IFdBc/20);
                    Lo2IF_q = loQ / 10^(Lo2IFdBc/20);
                    if 0
                        pwrdB(Lo2IF_i) - pwrdB(loI)
                    end

                    switch flag_MixType
                        case {'IQ', 'IMR'}
                            xI = xI + Lo2IF_i;
                            xQ = xQ + Lo2IF_q;
                        case {'Real', 'SSB'}
                            xIF = xIF + Lo2IF;
                    end
                end
            end

            % iq imbalance pre-compesnation
            if  obj.flag_IQimb_comp && strcmpi(obj.flag_UDconvert, 'U')
                xI_comp = 10^(-0.5*obj.imb_AmpDB_comp/20)*xI.*exp(-1i*0.5*obj.imb_PhsDeg_comp/180*pi);
                xQ_comp = 10^(+0.5*obj.imb_AmpDB_comp/20)*xQ.*exp(+1i*0.5*obj.imb_PhsDeg_comp/180*pi);

                xI = xI_comp;
                xQ = xQ_comp;
            end
            % iq imbalance
            if flag_IQimbalance && strcmpi(flag_MixType, 'IQ')
                methodIQimb = 'kc';
                switch methodIQimb
                    case 'matlab'
                        xIQ_imb2 = iqimbal(xI+1i*xQ,obj.imb_AmpDB,obj.imb_PhsDeg);
                        xI = xIQ_imb2(:,1);
                        xQ = xIQ_imb2(:,2);
                    case 'kc'
                        xIQ_imb = obj.iqImbalance([xI, xQ], obj.imb_AmpDB, obj.imb_PhsDeg);
                        xI = xIQ_imb(:,1);
                        xQ = xIQ_imb(:,2);
                end
                if 0
                    [IMB_MagDB_est, IMB_PhsDeg_est] = IQ_IMB_cor(xI,  xQ, 'TXQEC')
                    [IMB_MagDB_est2, IMB_PhsDeg_est2] = IQ_IMB_cor(xI_comp,  xQ_comp, 'TXQEC')
                    PLOT_FFT_dB(xI,fs,[],[],0510,'xI');
                    PLOT_FFT_dB(xIQ_imb(:,1),fs,[],[],0510,'xImbI');
                    PLOT_FFT_dB(xQ,fs,[],[],0510*10,'xQ');
                    PLOT_FFT_dB(xIQ_imb(:,2),fs,[],[],0510*10,'xImbQ');

                end
            end

            % lo leakage to RF
            if flag_LoLeakage
                Lo2RFdBc = abs(obj.Lo2RF_dBc);
                if Lo2RFdBc~=0
                    Lo2RF = lo / 10^(Lo2RFdBc/20);
                    Lo2RF_i = loI / 10^(Lo2RFdBc/20);
                    Lo2RF_q = loQ / 10^(Lo2RFdBc/20);
                    if 0
                        pwrdB(Lo2RF_i) - pwrdB(loI)
                    end

                    if flag_IQimbalance && strcmpi(flag_MixType, 'IQ')
                        [Lo2RF_imb] = obj.iqImbalance([Lo2RF_i, Lo2RF_q], obj.imb_AmpDB, obj.imb_PhsDeg);
                        Lo2RF_i= Lo2RF_imb(:,1);
                        Lo2RF_q = Lo2RF_imb(:,2);
                    end
                else
                    Lo2RF = 0;
                    Lo2RF_i = 0;
                    Lo2RF_q = 0;
                end
            else
                Lo2RF = 0;
                Lo2RF_i = 0;
                Lo2RF_q = 0;
            end

            % iq level offset pre-compensation, 2023-05-11
            if strcmpi(obj.flag_UDconvert, 'U') && obj.flag_IQLvlOffset_comp
                error('To compensate the IQ level offset, the source should be cw signal !')
            end

            % iq level offset
            if flag_IQLvlOffset && strcmpi(flag_MixType, 'IQ')
                LvlOffsetDb = obj.imb_LvlOffsetDB;
                if numel(LvlOffsetDb) == 2
                    LvlOffset_i = 10^(LvlOffsetDb(1)/20);
                    LvlOffset_q = 10^(LvlOffsetDb(2)/20);
                else
                    LvlOffset_i = 10^(+LvlOffsetDb/2/20);
                    LvlOffset_q = 10^(-LvlOffsetDb/2/20);
                end
                xI = xI + (LvlOffset_i-1);
                xQ = xQ + (LvlOffset_q-1);
                if 0
                    pwrdB(xI) - pwrdB(xQ)
                end
            end

            % mixing
            switch flag_MixType
                case 'IQ'
                    switch flag_UDconvert
                        case 'U'
                            yI = xI.*loI + Lo2RF_i;
                            yQ = xQ.*loQ + Lo2RF_q;
                        case 'D'
                            if 0
                                if strcmpi(obj.typeInputSignal,'complx')
                                    error('mixer.mix.flag_UDconvert.typeInputSignal: Not suppor complex signal for downconversion!')
                                end
                            end
                            yI = xIF.*loI + Lo2RF_i;
                            yQ = xIF.*loQ + Lo2RF_q;
                        otherwise
                            error('mixer.mix.flag_UDconvert: check the input of flag_UDconvert ?!')
                    end
                    if exist('phaseDiffCal.m') ~= 2 % check fucntion (phaseDiffCal)
                        error('mixer.mix: check the funciton phaseDiffCal.m exist ?!')
                    end

                    typeWaveform = 'carrier';
                    switch typeWaveform
                        case 'cw'
                            phsDiff_x = phaseDiffCal(xI, xQ);
                            phsDiff_lo = phaseDiffCal(loI, loQ);
                            if phsDiff_lo == 0
                                error('mixer.mix: check the output phsDiff_lo ?!')
                            end

                            a = fix(phsDiff_x/phsDiff_lo);
                            if a==0
                                a = fix(phsDiff_lo/phsDiff_x);
                            end
                        otherwise
                            a = 1;
                    end

                    if strcmpi(obj.flag_UDconvert, 'U')
                        switch obj.flag_UDconvert_FreqSelect
                            case 'H'
                                B = -1 * a;
                                y = yI + B*yQ;
                            case 'L'
                                B = +1 * a;
                                y = yI + B*yQ;
                        end
                    elseif strcmpi(obj.flag_UDconvert, 'D')
                        switch obj.flag_UDconvert_FreqSelect
                            case 'H'
                                B = -1 * a;
                            case 'L'
                                B = +1 * a;
                        end
                        y = yI + B*1i*yQ;
                    end

                case {'Real','REAL','real'}
                    y = xIF.*lo(:,1) + Lo2RF;
                    yI = [];
                    yQ = [];
            end

            if strcmpi(obj.flag_UDconvert, 'D')
                if ~isempty(obj.IFfilterBW)
                    flag_IFfilter = 1;
                    if isempty(obj.IFfilterTyp)
                        obj.IFfilterTyp = 'lpf';
                    end
                    [bIFfir, yIFfir] = firGen(y, fs, obj.IFfilterBW, 1e6,  obj.IFfilterTyp, 'eqrip', [0.001 50], []);
                    y = yIFfir;
                    yI = conv(bIFfir, yI, 'same');
                    yQ = conv(bIFfir, yQ, 'same');
                else
                    flag_IFfilter = 0;
                end

                if obj.fNCO ~= 0
                    if 0
                        if ~flag_IFfilter % 2023-05-10, add IF filter before ADC process(Decimation/NCO/...)
                            if isempty(obj.IFfilterBW)
                                obj.IFfilterBW = 0.5*(obj.fs/2 - 10e6)*[-1 1];
                            end
                            if isempty(obj.IFfilterTyp)
                                obj.IFfilterTyp = 'lpf';
                            end
                            [bIFfir, yIFfir] = firGen(y, fs, obj.IFfilterBW, 1e6,  obj.IFfilterTyp, 'eqrip', [0.1 50], []);
                            y = yIFfir;
                        end
                    end
                    [yNCO, yIQ] = nco(y, obj.fNCO, fs, [], []);
                    y = yNCO;
                    yI = yIQ(:,1);
                    yQ = yIQ(:,2);
                end
            end

            % power cal.
            obj.pwrdBm = pwrdB(y,[],[],[],'dBm');
            if 0
                pwrdBmY_Cal = 10*log10(pwrdB(xIF,[],[],[],'W')*pwrdB(lo,[],[],[],'W')*1000)
                pwrdBmLo = pwrdB(lo,[],[],[],'dBm')
                pwrdBmXif = pwrdB(xIF,[],[],[],'dBm')
                pwrdBmY = pwrdB(y,[],[],[],'dBm')
            end

            % export and update object
            if flag_row
                obj.y = y.';
                obj.yI = yI.';
                obj.yQ = yQ.';
            else
                obj.y = y;
                obj.yI = yI;
                obj.yQ = yQ;
            end

            % plot
            if ~isempty(fnum)
                obj.plot(y, fs, 'dBm', 'kHz', 0, 'mix product', fnum);
                if 0
                    obj.plot(y, fs, 'dBm', 'kHz', 0, 'nan', [obj.fnum,3,1,1]);
                    obj.plot(y, fs, 'dBm', 'kHz', 0, 'lo leakage', [obj.fnum,3,1,2]);
                    obj.plot(y, fs, 'dBm', 'kHz', 0, 'imbalance', [obj.fnum,3,1,3]);
                    PLOT_FFT_dB(y,fs,[],[],123);
                end
            end

        end

        function y = iqImbalance(obj, x, imb_AmpDB, imb_PhsDeg)
            if exist('iqImbal.m') ~= 2 % check fucntion (iqImbal)
                error('mixer.iqImbalance: check the funciton iqImbal.m exist ?!')
            end
       
            if ~exist('imb_AmpDB', 'var') || isempty(imb_AmpDB)
                imb_AmpDB = obj.imb_AmpDB;
            end

            if ~exist('imb_PhsDeg', 'var') || isempty(imb_PhsDeg)
                imb_PhsDeg = obj.imb_PhsDeg;
            end

            % imbalance
            y = iqImbal(x, imb_AmpDB, imb_PhsDeg);

        end

        function plot(obj, x, fs, pwrUnit, freqUnit, freqLog, dispLegend, fnum)
            if exist('PLOT_FFT_dB.m') ~= 2 % check fucntion (plot)
                error('mixer.mix: check the funciton PLOT_FFT_dB.m exist ?!')
            end

            fs = obj.fs;
            if exist('fnum','var') && ~isempty(fnum) && ~isempty(x) && ~isempty(fs) 
                switch obj.flag_UDconvert_FreqSelect
                    case {'Zero'}
                        PLOT_FFT_dB(x, fs, [], pwrUnit, fnum, dispLegend, [], [], [], freqLog, freqUnit, 1, []);
                    otherwise
                        PLOT_FFT_dB(x, fs, [], pwrUnit, fnum, dispLegend, [], [], [], freqLog, freqUnit, 3, []);
                end
            end
        end
    end


end