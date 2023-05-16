function evm = emvNRdl_hRawEVM(meaSymbols, refSymbols)
%   hRawEVM Error vector magnitude calculation
%   Copyright 2020 The MathWorks, Inc.

%   EV   - The normalized error vector
%   EV=(X-R)/sqrt(mean(abs(R.^2))). This allows for peak and RMS EVM
%   calculation for pre-existing normalized error vectors.

if (nargin == 2)
    if iscell(meaSymbols)
        m = meaSymbols{:};
        r = refSymbols{:};
    else
        m = meaSymbols;
        r = refSymbols;
    end
    errorVector = m-r;
    p = sqrt(mean(abs(r(:).^2)));
    if (p == 0)
        p = 1;
    end
    evNorm = errorVector/p;
else
    
    if ~isempty(meaSymbols)
        ev = meaSymbols{1};
    else
        ev = [];
    end
    evNorm = ev;
end

evm.EV = evNorm;
evNormMag = abs(evNorm(:));
evm.EVNormMag = evNormMag;
evm.Peak = max(evNormMag);
evm.RMS = sqrt(mean(evNormMag.^2));

end
