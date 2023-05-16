function W = getEVMWindow(carrier,frequencyRange,channelBandwidth,nFFT)
%   W = getEVMWindow(CARRIER,FREQUENCYRANGE,CHANNELBANDWIDTH,NFFT) is the
%   error vector magnitude window length, as mentioned in TS 38.104, Section
%   B.5.2 < FR1 /FR2 >. W is defined for a given combination of subcarrier
%   spacing, channel bandwidth/fft length, frequency range and CP type. W
%   is subsequently used as an intermediate value to decide the CP Fraction
%   for OFDM demodulation

    scsFR1 = [15 30 60];
    scsFR2 = [60 120];
    % BW MHz        5  10 15 20  25  30  40  50  60  70  80  90  100
    nfftFR1 = [256 384 512 768 1024 1536 2048 3072  4096];
    WsFR1   = [NaN NaN  14 NaN   28   44   58  108   144;      % NormalCp, 15kHz
                 8 NaN  14  22   28   54   72  130   172;      % NormalCp, 30kHz
                 8  11  14  26   36   64   86  NaN   NaN;      % NormalCp, 60kHz
                54  80 106 164  220  340  454  NaN   NaN];     % ExtendedCp, 60kHz

    % BW MHz        50 100 200 400
    nfftFR2 = [512 1024 2048 4096];
    WsFR2        = [NaN  36  72  144;                           % NormalCP, 60kHz
                    18   36  72  144;                          % NormalCP, 120kHz
                   NaN  220  440 880];                         % ExtendedCP, 60kHz
    W = []; %#ok<NASGU>
    if (strcmpi(frequencyRange,'FR1'))
        rowIdx = find(carrier.SubcarrierSpacing == scsFR1) + double(strcmpi(carrier.CyclicPrefix,'extended'));
        W = WsFR1(rowIdx,nFFT == nfftFR1);
        if channelBandwidth == 25
            if nFFT == 512 && carrier.SubcarrierSpacing == 60
                if strcmpi(carrier.CyclicPrefix,'extended')
                    W = 110;
                else
                    W = 18;
                end
            elseif nFFT == 1024 && carrier.SubcarrierSpacing == 30
                W = 36;
            elseif nFFT == 2048 && carrier.SubcarrierSpacing == 15
                W = 72;
            end
        end
    else
        rowIdx = find(carrier.SubcarrierSpacing == scsFR2) + double(strcmpi(carrier.CyclicPrefix,'extended'))*2;
        W = WsFR2(rowIdx,nFFT == nfftFR2);
    end
    % Filter out invalid combinations
    if isnan(W) || isempty(W)
        error('Invalid FFT/SCS/BW combination');
    end
end