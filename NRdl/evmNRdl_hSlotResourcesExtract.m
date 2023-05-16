function [pdschIndices,pdschSymbols,dmrsIndices,dmrsSymbols,ptrsIndices,ptrsSymbols] = hSlotResourcesExtract(InfoResourcePDSCH, nSlot)
%   hSlotResources Slot resources extraction
%   Copyright 2019-2021 The MathWorks, Inc.

if 1 % Initialization
    NconfigPDSCH = length(InfoResourcePDSCH);
end

pdschIndices = [];
pdschSymbols = [];
dmrsIndices = [];
dmrsSymbols = [];
ptrsIndices = [];
ptrsSymbols = [];

for n = 1:NconfigPDSCH
    slotRange = [InfoResourcePDSCH(n).Resources.NSlot];
    [~,slotIdx] = ismember(nSlot,slotRange);
    if slotIdx % If the PDSCH is present in this slot, get the indices and symbols

        pdschIndices = [pdschIndices; InfoResourcePDSCH(n).Resources(slotIdx).ChannelIndices]; %#ok<*AGROW>
        pdschSymbols = [pdschSymbols; InfoResourcePDSCH(n).Resources(slotIdx).ChannelSymbols];
        dmrsIndices = [dmrsIndices; InfoResourcePDSCH(n).Resources(slotIdx).DMRSIndices];
        dmrsSymbols = [dmrsSymbols; InfoResourcePDSCH(n).Resources(slotIdx).DMRSSymbols];
        ptrsIndices = [ptrsIndices; InfoResourcePDSCH(n).Resources(slotIdx).PTRSIndices];
        ptrsSymbols = [ptrsSymbols; InfoResourcePDSCH(n).Resources(slotIdx).PTRSSymbols];
    end
end
end
