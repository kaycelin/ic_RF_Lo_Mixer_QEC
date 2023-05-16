%% 2022-11-16, draft
%% 2023-05-05, export yI and yQ for RxQEC

classdef Lo < handle
    properties
        % Lo sine waveform
        fs = []
        T = 1e-3
        fLo = []
        Nsamps = []
        flag_LoQuad = 'IQ' % 'IQ': [-45, 45] (degree), 'Real': 0 degree 

        % Lo impairment
        flag_LoQuad_method = 'LoIQ->LoPN' % 'LoIQ->LoPN'/'LoPN->LoIQ'
        flag_LoPhaseNoise = 0
        PwrdBm = 10
        phsNzFreqOffsetHz = [ 100e3, 200e3, 400e3, 600e3, 800e3,...
            1.2e6, 1.8e6, 6e6, 10e6 ]; % Offset From Carrier
        phsNzLvldBc = [ -107, -115, -123, -128, -131,...
            -136, -140, -151, -156 ]*1.0 + -0*10; % Phase Noise power

        % Lo waveform output
        y
        yI % 2023-05-05
        yQ % 2023-05-05

        % Plot
        fnum
    end

    properties (Dependent)

    end

    methods (Access = public)
        function [obj, obj_struct] = Lo(varargin)
            % Initialization and set parameters
            if iscell(varargin) && size(varargin,2)>1
                varargin_tmp = varargin;
            else
                varargin_tmp = varargin{:};
            end
            obj.setParameters(varargin_tmp);

            % generate Lo
            obj.LoGen();

            % export structure
            obj_struct = struct(obj);
        end
    end

    properties (Access = private)
        % flag_setParameters = 0, reset for setting
        flag_setParameters = 0
        default
        Nstages = 1
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
                error('Lo.setParameters: type of params is STRUCT or CELL!')
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
                    loParams = InputTmp{n,i};
                elseif isstruct(InputTmp)
                    loParams = InputTmp;
                else
                    error('Lo.setParameters: type of params type is STRUCT or CELL!')
                end

                if Nvar>i-1 && ~isempty(loParams)
                    if isstruct(loParams)
                        % lo parameters
                        if isfield(loParams,'fs')&&~isempty(loParams.fs)
                            obj.fs = loParams.fs;
                        end
                        if isfield(loParams,'T')&&~isempty(loParams.T)
                            obj.T = loParams.T;
                        end
                        if isfield(loParams,'fLo')&&~isempty(loParams.fLo)
                            obj.fLo = loParams.fLo;
                        end
                        if isfield(loParams,'Nsamps')&&~isempty(loParams.Nsamps)
                            obj.Nsamps = loParams.Nsamps;
                        end
                        if isfield(loParams,'flag_LoQuad')&&~isempty(loParams.flag_LoQuad)
                            obj.flag_LoQuad = loParams.flag_LoQuad;
                        end
                        if isfield(loParams,'flag_LoQuad_method')&&~isempty(loParams.flag_LoQuad_method)
                            obj.flag_LoQuad_method = loParams.flag_LoQuad_method;
                        end
                        if isfield(loParams,'flag_LoPhaseNoise')&&~isempty(loParams.flag_LoPhaseNoise)
                            obj.flag_LoPhaseNoise = loParams.flag_LoPhaseNoise;
                        end
                        if isfield(loParams,'PwrdBm')&&~isempty(loParams.PwrdBm)
                            obj.PwrdBm = loParams.PwrdBm;
                        end
                        if isfield(loParams,'phsNzFreqOffsetHz')&&~isempty(loParams.phsNzFreqOffsetHz)
                            obj.phsNzFreqOffsetHz = loParams.phsNzFreqOffsetHz;
                        end
                        if isfield(loParams,'phsNzLvldBc')&&~isempty(loParams.phsNzLvldBc)
                            obj.phsNzLvldBc = loParams.phsNzLvldBc;
                        end
                        if isfield(loParams,'fnum')&&~isempty(loParams.fnum)
                            obj.fnum = loParams.fnum;
                        end
                    else
                        error('Lo.setParameters: type of loParams type should be STRUCT!')
                    end
                end
                i = i+1;
            end
        end
    end

    methods
        function [y, pn] = LoGen(obj, fLo, PwrdBm, flag_LoQuad_method) % generate Lo
            % input parameters updated
            if ~exist('fLo', 'var') || isempty(fLo)
                fLo = obj.fLo;
            end
            if ~exist('PwrdBm', 'var') || isempty(PwrdBm)
                PwrdBm = obj.PwrdBm;
            end
            if ~exist('flag_LoQuad_method', 'var') || isempty(flag_LoQuad_method)
                flag_LoQuad_method = obj.flag_LoQuad_method;
            end
            if isempty(obj.fs)
                error('Lo.LoGen: check the sampling rate parameters ?!')
            end
            if isempty(obj.Nsamps)&&~isempty(obj.T)
                obj.Nsamps = obj.fs*obj.T;
            elseif isempty(obj.T)
                error('Lo.LoGen: check the Nsamps or T parameters ?!')
            end
            if isempty(fLo)
                error('Lo.LoGen: check the fLo parameters ?!')
            else
                t = [0:obj.Nsamps-1]/obj.fs; % column
                LoSignal_i = sqrt(2*0.001)*10^(PwrdBm/20)*sin(2*pi*fLo*t(:)); % generate real sinwave with 0dBm
                if 0
                    LoSignal_q = sqrt(2*0.001)*10^(PwrdBm/20)*cos(2*pi*fLo*t(:)); % generate real sinwave with 0dBm
                else
                    LoSignal_q = sqrt(2*0.001)*10^(PwrdBm/20)*sin(2*pi*fLo*t(:)+pi/2); % generate real sinwave with 0dBm
                end
            end

            switch flag_LoQuad_method
                case 'LoPN->LoIQ' % PhaseNoise apply then Quadurature applied
                    if obj.flag_LoPhaseNoise
                        [LoSignalPN, pn] = obj.LoPhsNoise(LoSignal, obj.fs, obj.phsNzFreqOffsetHz, obj.phsNzLvldBc, []);
                    else
                        LoSignalPN = LoSignal;
                        pn = [];
                    end

                    if exist('phaseShiftTime.m') ~= 2
                        error('Lo.LoGen: check the funciton phaseShiftTime.m exist ?!')
                    end
                    LoSignalPN_i = phaseShiftTime(LoSignalPN, 45, obj.fs, fLo); % gnerate i
                    LoSignalPN_q = phaseShiftTime(LoSignalPN, -45, obj.fs, fLo); % generate Q
                case 'LoIQ->LoPN' % Quadurature apply then PhaseNoise applied
                    LoSignal_i = cwGen(obj.fs, obj.Nsamps, PwrdBm, fLo, +45, []);
                    LoSignal_q = cwGen(obj.fs, obj.Nsamps, PwrdBm, fLo, -45, []);
                    if obj.flag_LoPhaseNoise
                        [LoSignalPN_i, pn_i] = obj.LoPhsNoise(LoSignal_i, obj.fs, obj.phsNzFreqOffsetHz, obj.phsNzLvldBc);
                        [LoSignalPN_q, pn_q] = obj.LoPhsNoise(LoSignal_q, obj.fs, obj.phsNzFreqOffsetHz, obj.phsNzLvldBc);
                        pn = [pn_i, pn_q];
                    else
                        LoSignalPN_i = LoSignal_i;
                        LoSignalPN_q = LoSignal_q;
                        pn = [];
                    end
                otherwise
                    error('Lo.LoGen: check the flag_LoQuad_method type ?!')
            end

            % export
            switch obj.flag_LoQuad
                case {'IQ','iq'}
                    y = [LoSignalPN_i, LoSignalPN_q];
                case {'Real','real'}
                    y = LoSignalPN_i;
            end

            % update object
            obj.y = y;

            % plot
            if ~isempty(obj.fnum)
                if 0
                    PLOT_FFT_dB(y(:,1), obj.fs, [], 'dBm', fnum_y, 'xi w/ phase noise', [], [], [], 1, 'kHz');
                    PLOT_FFT_dB(pn, obj.fs, [], 'dB', fnum_pn, 'phase noise', [], [], [], 1, 'kHz');
                else
                    freqMarkerPN = 0+[0, obj.phsNzFreqOffsetHz];
                    freqMarkerY = fLo+[0, obj.phsNzFreqOffsetHz];
                    freqScale = 1e3;
                    switch obj.flag_LoQuad
                        case {'IQ','iq'}
                            dispTitleLo = {'Lo I','Lo Q'};
                            dispTitlePN = {'PN I','PN Q'};
                            if obj.flag_LoPhaseNoise
                                fNUM = [obj.fnum, 2,2,0];
                            else
                                fNUM = [obj.fnum, 1,2,0];
                            end
                        case {'Real','real'}
                            dispTitleLo = {'Lo'};
                            dispTitlePN = {'PN'};
                            if obj.flag_LoPhaseNoise
                                fNUM = [obj.fnum, 1,2,0];
                            else
                                fNUM = [obj.fnum];
                            end
                    end
                    for k=1:size(obj.y, 2)
                        fNUM(end) = fNUM(end)+1;
                        obj.plot(y(:,k), obj.fs, 'dBm', freqScale, 0, {'xi w/ phase noise', dispTitleLo{k}}, {freqMarkerY}, fNUM);
                        xlim((fLo+[-1 1]*max(obj.phsNzFreqOffsetHz))/freqScale)

                        if obj.flag_LoPhaseNoise
                            fNUM(end) = fNUM(end)+1;
                            obj.plot(pn(:,k), obj.fs, 'dBc', freqScale, 1, {'phase noise', dispTitlePN{k}}, {freqMarkerPN}, fNUM);
                        end
                    end
                end
            end

        end

        function [y, pn] = LoPhsNoise(obj, x, fs, phsNzFreqOffsetHz, phsNzLvldBc, fnum)
            if ~exist('fs', 'var') || isempty(fs)
                fs = obj.fs;
            end

            if ~exist('phsNzFreqOffsetHz', 'var') || isempty(phsNzFreqOffsetHz)
                phsNzFreqOffsetHz = obj.phsNzFreqOffsetHz;
            end

            if ~exist('phsNzLvldBc', 'var') || isempty(phsNzLvldBc)
                phsNzLvldBc = obj.phsNzLvldBc;
            end

            if ~exist('fnum', 'var') || isempty(fnum)
                fnum = [];
            end

            Nsamps = numel(x);
            RBW = fs/Nsamps;
            
            methodPhsNZ = 'test';
            phsNzLvldBcHz = phsNzLvldBc - 10*log10(RBW);
            switch methodPhsNZ
                case 'matlab'
                    pnoise = comm.PhaseNoise('Level',phsNzLvldBcHz, ...
                        'FrequencyOffset',phsNzFreqOffsetHz, ...
                        'SampleRate',fs);
                    y = pnoise(x);
                    pn = [];
                case 'test'
                    [y, pn] = LoPhaseNoiseMatlab(x, fs, phsNzLvldBcHz, phsNzFreqOffsetHz);
                case 'test2'
                    error('Lo.LoPhsNoise: methodPhsNZ:test2, to be updated !')
                    [y, pn] = LoPhaseNoiseAlex(x, fs, phsNzFreqOffsetHz, phsNzLVLdBcHz, [], obj.fLo, Nsamps, obj.PwrdBm);
            end

            if ~isempty(fnum)
                if 1
                    scaleFreq = 1e3;
                    freqMarker = obj.fLo+[0, phsNzFreqOffsetHz];
                    mkPwrs = obj.plot(y, fs, 'dBm', scaleFreq, 1, 'xi w/ phase noise', freqMarker, fnum);
                    xlim([obj.fLo+[-1 1]*max(phsNzFreqOffsetHz)]/scaleFreq)
                    if 0
                        PwrdBcPerdBc = mkPwrs(2:end) - mkPwrs(1)
                        DiffPhNz = PwrdBcPerdBc - phsNzLvldBc
                        pwrdB(y,[],[],[],'dBm')
                    end
                end
            end
        end

        function mkPwrs = plot(obj, x, fs, pwrUnit, freqScale, freqLog, dispLegend, freqMarker, fnum)
            fs = obj.fs;
            if exist('fnum','var') && ~isempty(fnum) && ~isempty(x)
                [~,~,~,mkPwrs] = PLOT_FFT_dB(x, fs, [], 'dBC', fnum, dispLegend, [], [], [], freqLog, freqScale, 3, freqMarker);
            end
        end
    end
end

