function evmNRdownlink(dlRefwaveConfig, dlRefwaveInfo, dlResourceInfo, rxWaveform)
%   Copyright The MathWorks, Inc.
%   hNRDownlinkEVM

if nargin~=4
    error('evmNRdownlink: No enough input')
end
if ~exist('dlRefwaveConfig','var')||isempty(dlRefwaveConfig)
    error('evmNRdownlink: dlRefwaveConfig is empty')
elseif ~strcmpi(class(dlRefwaveConfig), 'nrDLCarrierConfig')
    error('evmNRdownlink: dlRefwaveConfig is not class of nrDLCarrierConfig')
elseif ~exist('dlRefwaveInfo','var')||isempty(dlRefwaveInfo)
    error('evmNRdownlink: dlRefwaveInfo is empty')
elseif ~isfield(dlRefwaveInfo, 'Info')
    error('evmNRdownlink: dlRefwaveInfo is not included of Info')
elseif ~exist('dlResourceInfo','var')||isempty(dlResourceInfo)
    error('evmNRdownlink: dlResourceInfo is empty')
elseif ~isfield(dlResourceInfo, 'WaveformResources')
    error('evmNRdownlink: dlResourceInfo is not included of WaveformResources')
elseif ~exist('rxWaveform','var')||isempty(rxWaveform)
    error('evmNRdownlink: waveform is empty')
end

if 1 % Initialization
    evm3GPP = 0;
    Nframes = 1;
    pdschCfgLen = 1;
    pdcchCfgLen = 0;
    displayEVM = 1;
    plotEVM = 1;
end

if 1 % Demodulation
    if 1 % Initializaion
        if 1 % Parameters
            Config = dlRefwaveConfig;
            Info = dlRefwaveInfo.Info;
            InfoResourcePDSCH = dlResourceInfo.WaveformResources.PDSCH;
        end

        if 1 % Set PDSCH configuration
            if strcmpi(class(Config.PDSCH{1}), 'nrWavegenPDSCHConfig')
                pdschCell = [];
                for idx = 1:length(InfoResourcePDSCH)
                    pdschCell{idx} = nrPDSCHConfig;
                    %                     pdsch = InfoResourcePDSCH(idx).PDSCH;
                    %                     pdsch = Config.PDSCH{1};
                    pdschCell{idx}.PRBSet = Config.PDSCH{idx}.PRBSet;
                    pdschCell{idx}.SymbolAllocation     = Config.PDSCH{idx}.SymbolAllocation;
                    pdschCell{idx}.Modulation           = Config.PDSCH{idx}.Modulation;
                    pdschCell{idx}.NumLayers            = Config.PDSCH{idx}.NumLayers;
                    pdschCell{idx}.MappingType          = Config.PDSCH{idx}.MappingType;
                    pdschCell{idx}.RNTI                 = Config.PDSCH{idx}.RNTI;
                    pdschCell{idx}.NID                  = Config.PDSCH{idx}.NID;
                    pdschCell{idx}.VRBToPRBInterleaving = Config.PDSCH{idx}.VRBToPRBInterleaving;
                    pdschCell{idx}.VRBBundleSize        = Config.PDSCH{idx}.VRBBundleSize;
                    for idx2 = 1:length(Config.PDSCH{idx}.ReservedPRB)
                        pdschObj{idx}.ReservedPRB{idx2} = nrPDSCHReservedConfig;
                        pdschObj{idx}.ReservedPRB{idx2}.PRBSet    = Config.PDSCH{idx}.ReservedPRB{idx2}.PRBSet;
                        pdschObj{idx}.ReservedPRB{idx2}.SymbolSet = Config.PDSCH{idx}.ReservedPRB{idx2}.SymbolSet;
                        pdschObj{idx}.ReservedPRB{idx2}.Period    = Config.PDSCH{idx}.ReservedPRB{idx2}.Period;
                    end
                    % Set DM-RS parameters
                    pdschCell{idx}.DMRS = Config.PDSCH{idx}.DMRS;

                    % Set PT-RS parameters
                    pdschCell{idx}.EnablePTRS = Config.PDSCH{idx}.EnablePTRS;
                    pdschCell{idx}.PTRS = Config.PDSCH{idx}.PTRS;
                    if Config.PDSCH{idx}.Enable && sum(Config.PDSCH{idx}.DMRS.CDMLengths)
                        cdmLengths = Config.PDSCH{idx}.DMRS.CDMLengths;
                    end
                end
            else
                error('evmNRdownlink: Config.PDSCH{1} is not class of nrWavegenPDSCHConfig')
            end
        end

        if 1 % Get parameters
            sampleRate = Info.SampleRate;
            initialNSlot = 0;
        end
    end

    for bwpIdx = 1:length(Config.BandwidthParts) % Store the received waveform for reuse in each BWP

        if 1 % Initialization
            if 1 % Set Carrier configuration, Carrier based on bwpIdx
                for carIdx = 1:length(Config.SCSCarriers)
                    Carrier = nrCarrierConfig;
                    Carrier.NCellID = carIdx;
                    Carrier.SubcarrierSpacing = Config.SCSCarriers{carIdx}.SubcarrierSpacing;
                    Carrier.CyclicPrefix = Config.BandwidthParts{bwpIdx}.CyclicPrefix;
                    Carrier.NSizeGrid = Config.SCSCarriers{carIdx}.NSizeGrid;
                    Carrier.NStartGrid = Config.SCSCarriers{carIdx}.NStartGrid;
                end
            end
            if 1
                NSymbsPerSlot = Carrier.SymbolsPerSlot;
                NSlots = Nframes * Carrier.SlotsPerFrame * Carrier.SlotsPerSubframe; % Num. of slots of carrier
                NSymbols = NSymbsPerSlot * NSlots; % Num. of OFDM symbols in ref. grid
                NSubcarriers = Carrier.NSizeGrid * 12; % Num. of subcarriers in ref. grid
                NLayers = Config.PDSCH{1}.NumLayers; % Num of layers in ref. grid
            end

        end

        if 1 % Generate a reference grid of length two frames for timing synchronization ?, hReferenceGrid

            if 1 % Initialization
                if ~exist('Config','var')||isempty(Config)
                    error('evmNRdownlink: hReferenceGrid miss input parameter of Config')
                else
                    NSizeBWP = Config.BandwidthParts{bwpIdx}.NSizeBWP;
                    NStartGridBWP = Config.BandwidthParts{bwpIdx}.NStartBWP-Carrier.NStartGrid;
                    indBWP = 12*NStartGridBWP+1:12*(NStartGridBWP+NSizeBWP);
                    NREsPerSlot = NSizeBWP*12*NSymbsPerSlot; % ?

                    bwpGridDMRS = zeros(NSizeBWP*12,NSymbols,NLayers); %
                    refGrid = zeros(NSubcarriers,NSymbols,NLayers); % reference grid
                end

                if ~exist('InfoResourcePDSCH','var')||isempty(InfoResourcePDSCH)
                    error('evmNRdownlink: hReferenceGrid miss input parameter of InfoResourcePDSCH')
                elseif ~isfield(InfoResourcePDSCH, 'Resources')
                    error('evmNRdownlink: hReferenceGrid, InfoResourcePDSCH no field of Resources')
                else % slotRange
                    slotRange = [];
                    for pIdx = 1:length(InfoResourcePDSCH)
                        slotRange = [slotRange InfoResourcePDSCH(pIdx).Resources.NSlot];
                    end
                    slotRange = unique(slotRange);
                end
                if ~exist('Carrier','var')||isempty(Carrier)
                    error('evmNRdownlink: hReferenceGrid miss input parameter of Carrier')
                else
                    Carrier.NSlot = initialNSlot;
                end
            end

            if 1 % DMRS and RefGrid mapping
                for nSlot = Carrier.NSlot + (0:NSlots-1)
                    isDataSlot = ismember(nSlot, slotRange);
                    if isDataSlot
                        [~,~,dmrsIndices,dmrsSymbols,~,~] = evmNRdl_hSlotResourcesExtract(InfoResourcePDSCH, nSlot);
                    end

                    if ~isempty(dmrsIndices)
                        for n = 1:NLayers
                            dmrsIndices(:,n) = dmrsIndices(:,n) - NREsPerSlot*(n-1) + (NSymbols*NSizeBWP*12*(n-1));
                            pwr = 1; % Power adjustment is needed only for PUSCH case
                            bwpGridDMRS(dmrsIndices(:,n)+(nSlot-Carrier.NSlot)*NREsPerSlot) = dmrsSymbols(:,n)*pwr;
                            refGrid(indBWP,:,:) = bwpGridDMRS;
                        end
                    end
                end
            end
        end

        if 1 % Frequency shifting, taking into account 'k0' from BWP
            if ~exist('Info','var')||isempty(Info)
                error('evmNRdownlink: Freq. shifting miss input parameter of Info')
            else
                sampleRate = sampleRate;
                k0 = Info.k0;
            end
            if ~exist('Carrier','var')||isempty(Carrier)
                error('evmNRdownlink: Freq. shifting miss input parameter of Carrier')
            end

            t = (0:(numel(rxWaveform)-1)).'/sampleRate;
            k0Offset = k0 * Carrier.SubcarrierSpacing * 1e3;
            rxWaveformk0Shifted = rxWaveform.*exp(-1i*2*pi*k0Offset*t);
        end

        if 1 % Time synchronization of input waveform
            if ~exist('rxWaveformk0Shifted','var')||isempty(rxWaveformk0Shifted)
                error('evmNRdownlink: Time synchronization should be do after Frequency shifting')
            elseif ~exist('refGrid','var')||isempty(refGrid)
                error('evmNRdownlink: Time synchronization miss the refGrid to calculate timing offset')
            elseif ~exist('Carrier','var')||isempty(Carrier)
                error('evmNRdownlink: Time synchronization, nrTimingEstimate miss the Carrier')
            else
                offset = nrTimingEstimate(Carrier,rxWaveformk0Shifted,refGrid,'SampleRate',sampleRate);
                rxWaveform = rxWaveformk0Shifted(1+offset:end,:);
            end
        end

        if 1 % Demouldate waveform into a rx grid
            if 1 % Initialization
                if ~exist('Carrier','var')||isempty(Carrier)
                    error('evmNRdownlink: Demoudulation miss the parameters of Carrier')
                end
                if ~exist('Config','var')||isempty(Config)
                    error('evmNRdownlink: Demoudulation miss the parameters of Config')
                elseif ~any(ismember(properties(Config), 'FrequencyRange'))
                    error('evmNRdownlink: Demoudulation, Config miss the field of FrequencyRange')
                elseif ~any(ismember(properties(Config), 'ChannelBandwidth'))
                    error('evmNRdownlink: Demoudulation, Config miss the field of ChannelBandwidth')
                end
                if ~exist('Info','var')||isempty(Info)
                    error('evmNRdownlink: Demoudulation miss the parameters of Info')
                elseif ~isfield(Info, 'Nfft')
                    error('evmNRdownlink: Demoudulation, Info miss the field of Nfft')
                elseif ~isfield(Info, 'SampleRate')
                    error('evmNRdownlink: Demoudulation, Info miss the field of SampleRate')
                end
                if ~exist('indBWP','var')||isempty(indBWP)
                    error('evmNRdownlink: Demoudulation miss the parameters of indBWP for BWP dimesion mapping')
                end
            end
            % When evm3GPP is true, two EVM window locations and two CP fractions
            % are selected for 3GPP EVM for OFDM demodulation. If false, a single
            % EVM window location is used, which is centred in the middle of the CP
            if evm3GPP
                W = evmNRdl_getEVMWindow(Carrier,Config.FrequencyRange,Config.ChannelBandwidth,Info.Nfft);
                nEVM = 2;
                cpFraction = [0 ; W/cpLength];
            else
                nEVM = 1;
                cpFraction = 0.5;      % Use default value
            end

            % Demodulate the waveform
            Carrier.NSlot = initialNSlot;
            rxGridLow = nrOFDMDemodulate(Carrier, rxWaveform, 'CyclicPrefixFraction',cpFraction(1),'SampleRate',Info.SampleRate);
            if evm3GPP
                rxGridHigh = nrOFDMDemodulate(Carrier, rxWaveform, 'CyclicPrefixFraction',cpFraction(2),'SampleRate',Info.SampleRate);
            end

            % Resize rxGridLow based on BWP dimensions, Work only on the relevant BWP in the waveform to simplify indexing
            rxGridLow = rxGridLow(indBWP,:,:);
            if evm3GPP
                rxGridHigh = rxGridHigh(indBWP,:,:);
            end
        end

        if 1 % Generate a reference grid for channel estimation.
            if 1
                if ~exist('Carrier','var')||isempty(Carrier)
                    error('evmNRdownlink: refGrid_channel_estimation misss the parameters of Carrier to calculate NSlotsGrid ')
                else
                    NSymbsPerSlot = Carrier.SymbolsPerSlot;
                end
                if ~exist('Info','var')||isempty(Info)
                    error('evmNRdownlink: refGrid_channel_estimation misss the parameters of Info to calculate NFramesGrid')
                end
                if ~exist('indBWP','var')||isempty(indBWP)
                    error('evmNRdownlink: refGrid_channel_estimation miss the parameters of indBWP for BWP dimesion mapping')
                end
                if ~exist('rxGridLow','var')||isempty(rxGridLow)
                    error('evmNRdownlink: refGrid_channel_estimation miss the parameters of rxGridLow for refGridChEst transformation')
                end
                if 1 % Calculate number of slots and frames for the gives sampleRate
                    NSlotsGrid = floor(size(rxGridLow,2)/NSymbsPerSlot);
                    NFramesGrid = floor(NSlotsGrid/(10*Info.SlotsPerSubframe));
                end
            end

            % Generate a reference grid, refGrid, for slots corresponding to
            % the length of the input waveform. This grid contains only the
            % DM-RS and is primarily used for channel estimation.
            refGridChEst = refGrid(:,1:NSymbsPerSlot*NSlotsGrid);

            % Resize refGrid based on BWP dimensions
            refGridChEst = refGridChEst(indBWP,:,:);
        end

        if 1 % hChannelEstEVM
            if 1 % Initialization
                if ~exist('evm3GPP','var')||isempty(evm3GPP)
                    error('evmNRdownlink: hChannelEstEVM miss the parameters of evm3GPP')
                end
            end

            if evm3GPP
                if ~exist('NSlots10Ms','var')||isempty(NSlots10Ms)
                    if ~exist('NSlotsGrid','var')||isempty(NSlotsGrid)
                        error('evmNRdownlink: hChannelEstEVM, evm3GPP miss the parameters of NSlotsGrid for NSlots10Ms')
                    elseif ~exist('Carrier','var')||isempty(Carrier)
                        error('evmNRdownlink: hChannelEstEVM, evm3GPP miss the parameters of Carrier.SubcarrierSpacing for NSlots10Ms')
                    else
                        NSlots10Ms = min(NSlotsGrid,10*Carrier.SubcarrierSpacing/15);
                    end
                end
                if ~exist('rxGridLow','var')||isempty(rxGridLow)
                    error('evmNRdownlink: hChannelEstEVM, nrChannelEstimate miss the parameters of rxGridLow')
                elseif ~exist('refGridChEst','var')||isempty(refGridChEst)
                    error('evmNRdownlink: hChannelEstEVM, nrChannelEstimate miss the parameters of refGridChEst')
                elseif ~exist('pdschCell','var')||isempty(pdschCell)
                    error('evmNRdownlink: hChannelEstEVM, nrChannelEstimate miss the parameters of pdschCell')
                end

                NBlocks = ceil(NSlotsGrid/NSlots10Ms);
                % For downlink 3GPP channel estimation, the averaging needs to be
                % done in chunks of 10ms. For uplink channel estimation, averaging
                % is done on a per-lsot basis.
                if ~exist('rxGridLow','var')||isempty(rxGridLow)
                    dlFlag = 0;
                    NIterations = NSlotsGrid;
                else
                    dlFlag = 1;
                    NIterations = NBlocks;
                end
                for idx = 1:NIterations
                    if dlFlag
                        symIdx = (idx-1)*NSymbsPerSlot*NSlots10Ms+(1:(NSymbsPerSlot*NSlots10Ms));
                    else
                        symIdx = (idx-1)*NSymbsPerSlot+1:idx*NSymbsPerSlot;
                    end
                    % Use a smoothing filter in the frequency direction when dlFlag is true
                    if dlFlag
                        HestLowBlk = hChannelEstimateEVM3GPP(rxGridLow(:, symIdx, :),refGrid(:, symIdx,:),'movingAvgFilter','CDMLengths',cdmLengths);
                        HestHighBlk = hChannelEstimateEVM3GPP(rxGridHigh(:, symIdx, :),refGrid(:, symIdx,:),'movingAvgFilter','CDMLengths',cdmLengths);
                    else
                        HestLowBlk = hChannelEstimateEVM3GPP(rxGridLow(:, symIdx, :),refGrid(:, symIdx,:),'CDMLengths',cdmLengths);
                        HestHighBlk = hChannelEstimateEVM3GPP(rxGridHigh(:, symIdx, :),refGrid(:, symIdx,:),'CDMLengths',cdmLengths);
                    end
                    HestLow  = [HestLow HestLowBlk]; %#ok<*AGROW>
                    HestHigh = [HestHigh HestHighBlk];
                    estChannelGrid = HestHigh;
                end
            else
                estChannelGrid = [];
                % Compute channel estimates for each slot
                for slotIdx = 1:NSlotsGrid
                    symbsIdx = (slotIdx-1)*NSymbsPerSlot+1:slotIdx*NSymbsPerSlot;
                    symbsIdx(symbsIdx>size(rxGridLow, 2)) = [];
                    estChannelGridSlot = nrChannelEstimate(rxGridLow(:,symbsIdx,:),refGridChEst(:,symbsIdx,:),'CDMLengths',cdmLengths);
                    estChannelGrid = [estChannelGrid estChannelGridSlot];
                end
            end
        end

        if 1 % hDecodeSlots
            if 1 % Initialization
                % Slot allocation of the PDSCH configurations may overlap with each
                % other. Extract unique allocated slots
                slotRange = [];
                for pIdx = 1:length(InfoResourcePDSCH)
                    slotRange = [slotRange InfoResourcePDSCH(pIdx).Resources.NSlot];
                end
                slotRange = unique(slotRange);
                % In case of non-3GPP case, only a single EVM grid is used.
                % Compute the DL EVM for each active/valid DL slot, store the results
                % in a cell-array for later processing. Skip slots which are not DL.
                slotRange = slotRange(slotRange < initialNSlot+NSlotsGrid);
                slotRange = slotRange(slotRange >= initialNSlot);

                eqSym = cell(nEVM,1);                  % Equalized symbols for constellation plot, for each low/high EVM window location
                refSym = cell(nEVM,1);                 % Reference symbols for constellation plot, for each low/high EVM window location
            end

            for slotIdx=slotRange
                for e = 1:nEVM
                    if 1 % refSym, based on pdschSymbols_refSym
                        [pdschIndices_refSym,pdschSymbols_refSym,...
                            dmrsIndices_refSym,dmrsSymbols_refSym,...
                            ptrsIndices_refSym,ptrsSymbols_refSym] = evmNRdl_hSlotResourcesExtract(InfoResourcePDSCH, slotIdx);

                        if 1
                            % pdschEncodingOn = true to enables derate matching and LDPC decoding.
                            pdschEncodingOn = false;
                            % skipCtrlSymbols = true to skips symbols occupied by PDCCH when measuring PDSCH EVM.
                            skipCtrlSymbols = true;
                            decodePdsch = true;
                        end

                        % Do not include first two slot symbols for PDSCH EVM (used for
                        % control as specified in TS 38.141-1 table 4.9.2.2-2 (NR-TMs) /
                        % TS 38.101-1 table A.3.1-1 (FRCs))
                        % This step is not performed when
                        % * extrapolateHest is true, or
                        % * label is absent, or
                        % * coding is enabled
                        % Ensure that the symbol allocation does not include the first two
                        % OFDM symbols when the 'label' field is set. A non-empty 'label'
                        % is assumed to contain a string associated with Release 15
                        % NR-TMs/FRCs
                        if decodePdsch && skipCtrlSymbols  && ~pdschEncodingOn
                            carIdx = pdschIndices_refSym(:,1) <= (2*NSubcarriers);
                            pdschIndices_refSym(carIdx,:) = [];
                            pdschSymbols_refSym(carIdx,:) = [];
                        end
                    end

                    % Accumulate reference IQs, refSym
                    refSym{e} = [refSym{e}; pdschSymbols_refSym];
                end

                if 1 % eqSym, equalized symbols based on rxGridLow and estChannelGrid
                    % Get PDSCH resource elements from the received grid
                    [pdschRx,pdschHest] = nrExtractResources(pdschIndices_refSym,rxGridLow,estChannelGrid);

                    noiseEst = 0; % ZF based equalization, as per TS 38.104, Annex B.1 (FR1) / Annex C.1 (FR2)

                    % Equalization
                    [pdschEqGrid,csiLow] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

                    % Accumulate equalized IQs, eqSym
                    eqSym{e} = [eqSym{e}; pdschEqGrid];
                end

                if 1 % Map each reference-slot, equalized slot to its correct position in grid, layer, EVM-edge
                    if 1 % Initialization
                        decodePusch = false;
                        tmpGrid = zeros(NSubcarriers,NSymbsPerSlot);
                        refSlotGrid = zeros(size(estChannelGrid,1),NSlotsGrid*NSymbsPerSlot,size(estChannelGrid,3),nEVM);
                        eqSlotGrid = refSlotGrid;
                    end

                    for layerIdx = 1:NLayers
                        if layerIdx <= size(pdschIndices_refSym,2)
                            if decodePusch
                                ind = zeros(size(puschIndCurrentSlot,1),NLayers);
                                ind(:,layerIdx) = puschIndCurrentSlot(:,layerIdx) - NREsPerSlot*(layerIdx -1);
                            else
                                ind = zeros(size(pdschIndices_refSym,1),NLayers);
                                ind(:,layerIdx) = pdschIndices_refSym(:,layerIdx) - NREsPerSlot*(layerIdx -1);
                            end
                            tmpGrid(ind(:,layerIdx)) = pdschSymbols_refSym(:,layerIdx);
                            refSlotGrid(:,symbsIdx,layerIdx,nEVM) = tmpGrid;
                            tmpGrid(:) = 0;
                            tmpGrid(ind(:,layerIdx)) = pdschEqGrid(:,layerIdx);
                            eqSlotGrid(:,symbsIdx,layerIdx,nEVM) = tmpGrid;
                            tmpGrid(:) = 0;
                        end
                    end
                end
            end
        end

        if 1 % calculate evm
            if 0
                evmTmp = hEVM(Carrier,eqSlotGrid,refSlotGrid);
                evm = evmTmp.EVM;
            else % hEVM
                if 1 % Initialization
                    rxSlotGrid = eqSlotGrid;

                    % Locate allocated symbols in grid
                    symRange = find(sum(rxSlotGrid,[1 3 4]))-1;
                    slotRange = unique(floor(symRange/NSymbsPerSlot));
                    if iscolumn(slotRange)
                        slotRange = slotRange.';
                    end

                    if 0
                        evm = repmat(hRawEVM([]), nEVM , NSlotsGrid);
                    else
                        evm = repmat(emvNRdl_hRawEVM([]), nEVM , NSlotsGrid);
                    end
                end

                if 1 % Calculate raw evm
                    % Loop over each valid slot with these steps
                    % - Extract allocated REs, reference IQs per layer, per edge
                    % - Calculate raw error vector
                    % - Construct up to 2 EVM grids corresponding to the allocated slots
                    for slotIdx = slotRange
                        for e = 1:nEVM
                            % Initialization
                            refSymbols = reshape(refSlotGrid(:,symbsIdx,:,e),NSubcarriers*NSymbsPerSlot,NLayers);
                            rxSymbols = reshape(eqSlotGrid(:,symbsIdx,:,e),NSubcarriers*NSymbsPerSlot,NLayers);

                            ind  = find(rxSymbols);
                            ind = reshape(ind,length(ind)/NLayers,NLayers);
                            rxSymbols = rxSymbols(ind);
                            refSymbols = refSymbols(ind);
                            if 0
                                evmTmp = hRawEVM(rxSymbols,refSymbols);
                                evm(e,slotIdx+1) = evmTmp;
                            else
                                evmTmp = emvNRdl_hRawEVM(rxSymbols, refSymbols);
                                evm(e,slotIdx+1) = evmTmp;
                            end

                            % Initialization, Build low and high edge EVM grids across all slots
                            evmSlotEdge = zeros(2,NSubcarriers,NSymbsPerSlot);

                            evmSlotEdge(e,ind(:,1)) = mean(evm(e).EVNormMag*100,2);
                            symbsIdx2 = (slotIdx)*NSymbsPerSlot+1:(slotIdx+1)*NSymbsPerSlot;
                            evmGridEdge(e,:,symbsIdx2) =  evmSlotEdge(e,:,:);
                        end
                    end
                end

                if 1 % Summary evm
                    if 1
                        evmInfo.PDSCH = [];
                        evmInfo.PDCCH = [];
                    end

                    % The low or high edge timing is chosen for plotting
                    % automatically based on whichever has the largest RMS across
                    % all slots (Largest RMS chosen as mentioned in TS 38.101-1/2 , Annex F.6)
                    evmGrid = squeeze(evmGridEdge(1,:,:));
                    if evm3GPP
                        evmMaxLow = max([evm(1,:).RMS]);
                        evmMaxHigh = max([evm(2,:).RMS]);
                        if evmMaxHigh > evmMaxLow % replace by max evm value
                            evmGrid = squeeze(evmGridEdge(2,:,:));
                        end
                    end

                    % RMS and Peak EVM versus subcarrier
                    evmSubcarrierRMS = sqrt(sum(evmGrid.^2,2)./sum(evmGrid~=0,2));
                    evmSubcarrierPeak = max(evmGrid,[],2)./any(evmGrid,2);

                    % RMS and Peak EVM versus OFDM symbol
                    evmSymbolRMS = sqrt(sum(evmGrid.^2,1)./sum(evmGrid~=0,1)).';
                    evmSymbolPeak = (max(evmGrid,[],1)./any(evmGrid~=0,1)).';

                    % RMS and Peak EVM versus slot
                    evmSlotRMS = [];
                    evmSlotPeak = [];
                    for slotIdx = 1:NSlotsGrid
                        evmSlotRMSTmp = max(evm(1:end,slotIdx).RMS)*100;
                        evmSlotPeakTmp = max(evm(1:end,slotIdx).Peak)*100;
                        if isempty(evmSlotPeakTmp)
                            evmSlotPeakTmp = NaN;
                        end
                        evmSlotRMS= [evmSlotRMS;evmSlotRMSTmp]; %#ok<*AGROW>
                        evmSlotPeak = [evmSlotPeak;evmSlotPeakTmp];
                    end

                    % export
                    evmInfo.PDSCH.SubcarrierRMS = evmSubcarrierRMS;
                    evmInfo.PDSCH.SubcarrierPeak = evmSubcarrierPeak;
                    evmInfo.PDSCH.SymbolRMS = evmSymbolRMS;
                    evmInfo.PDSCH.SymbolPeak = evmSymbolPeak;
                    evmInfo.PDSCH.SlotRMS = evmSlotRMS;
                    evmInfo.PDSCH.SlotPeak = evmSlotPeak;
                    evmInfo.PDSCH.EVMGrid = evmGrid;
                    evmInfo.PDSCH.OverallEVM = evm;

                    % Below loop needs to run at least once
                    if (NFramesGrid == 0)
                        NFramesGrid = 1;
                    end

                    if 1
                        %      BandwidthPartID - Index of the bandwidth part
                        %      AverageFO       - Estimated frequency estimate
                        evmInfo.PDSCH.BandwidthPartID = Config.BandwidthParts{bwpIdx}.BandwidthPartID;
                        if 1 % averageFO
                            % When correctFineFO is enabled, estimate FO for each valid slot
                            % and store the averaged result in averageFO
                            averageFO = 0;
                        end
                        evmInfo.PDSCH.AverageFO = averageFO;
                    end
                end
            end
  
            if plotEVM
                % For each valid slot, update the plot (symbol,SC,slot,grid-wise)
                if pdschCfgLen
                    pIdx = find([evmInfo.PDSCH(:).BandwidthPartID] == Config.BandwidthParts{bwpIdx}.BandwidthPartID);
                    hEVMPlots(evmInfo.PDSCH(pIdx),eqSym{1,bwpIdx},refSym{1,bwpIdx},'PDSCH');
                end
                if pdcchCfgLen
                    pIdx = find([evmInfo.PDCCH(:).BandwidthPartID] == Config.BandwidthParts{bwpIdx}.BandwidthPartID);
                    hEVMPlots(evmInfo.PDCCH(pIdx),pdcchEqSym{1,bwpIdx},pdcchRefSym{1,bwpIdx},'PDCCH');
                end
            end

            if displayEVM
                if pdschCfgLen
                    pIdx = find([evmInfo.PDSCH(:).BandwidthPartID] == Config.BandwidthParts{bwpIdx}.BandwidthPartID);
                    fprintf('Averaged overall PDSCH RMS EVM: %0.3f%%\n', evmInfo.PDSCH(pIdx).OverallEVM.RMS*100);
                    disp("Overall PDSCH Peak EVM = " + string((evmInfo.PDSCH(pIdx).OverallEVM.Peak)*100) + "%");
                end
                if pdcchCfgLen
                    pIdx = find([evmInfo.PDCCH(:).BandwidthPartID] == Config.BandwidthParts{bwpIdx}.BandwidthPartID);
                    fprintf('Averaged overall PDCCH RMS EVM: %0.3f%%\n', evmInfo.PDCCH(pIdx).OverallEVM.RMS*100);
                    disp("Overall PDCCH Peak EVM = " + string((evmInfo.PDCCH(pIdx).OverallEVM.Peak)*100) + "%");
                end
            end
        end
    end
end
end
