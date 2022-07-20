function [dat,dF,dFOrg,opts] = baselineRemoveAndNoiseEstimation(datOrg,opts,evtSpatialMask,ff)
    
    % valid region
    [H,W,T] = size(datOrg);
    msk000 = var(datOrg,0,3)>1e-8;
    if exist('evtSpatialMask','var') && ~isempty(evtSpatialMask)
        evtSpatialMask = evtSpatialMask.*msk000;
    else
        evtSpatialMask = msk000;
    end
        
    % ignore the nan value
    nonvalid = isnan(datOrg);
    se = strel(true(ceil(2*opts.smoXY+1),ceil(2*opts.smoXY+1),5));
    nonvalid = imdilate(nonvalid,se);
    
    % smooth the data (memory-efficient version)
    dat = datOrg;
    if opts.smoXY>0
        for tt=1:size(dat,3)
            dat(:,:,tt) = imgaussfilt(dat(:,:,tt),opts.smoXY);
        end
    end
    
    %% linear estimation of F0
    opts.cut = min(opts.cut,T);
    [F0] = burst.baselineLinearEstimate(dat,opts.cut,opts.movAvgWin,nonvalid,ff);
    dFOrg = datOrg - F0;
    dF = dat - F0;
    
    % noise estimation
    stdMapGau = burst.noiseEstimate(dF,nonvalid,1,evtSpatialMask);
    stdMapGauOrg = burst.noiseEstimate(dFOrg,nonvalid,1,evtSpatialMask);
    if(opts.noiseEstimation == 1)
        stdMapGau(:) = median(stdMapGau(:));
        stdMapGauOrg(:) = median(stdMapGauOrg(:));
    end
    
    waitbar(0.4,ff,'Detecting active regions...');
    
    % normalization
    dF = dF./repmat(stdMapGau,1,1,T);
    dF = dF + opts.xBias;
    dFOrg = dFOrg./repmat(stdMapGauOrg,1,1,T);
    dFOrg = dFOrg + opts.xBias;
    opts.stdMap = stdMapGau;
end

