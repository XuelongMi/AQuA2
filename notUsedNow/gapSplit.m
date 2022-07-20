function [seLst,seLstInfoLabel] = gapSplit(subEvtLst,dat,dF,majorInfo,mergingInfo,ccRegions,opts)
    [H,W,T] = size(dF);
%     dff = dF.*repmat(stdMapGau,1,1,T)./(dat - dF.*repmat(stdMapGau,1,1,T));
    dff = dF./(dat./repmat(opts.stdMap,1,1,T) - dF);
    
    %% identify gap and consider cannot merge
    mergingInfo = gapIdentify2D(subEvtLst,dff,mergingInfo,majorInfo,opts);    
    [seLst,seLstInfoLabel] = se.mergingSEbyInfo(subEvtLst,majorInfo,mergingInfo,[H,W,T],ccRegions,opts);
end