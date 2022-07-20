function stdMapGau = noiseEstimate(xx,nonvalid,gap,evtSpatialMask)
T = size(xx,3);
xx(nonvalid) = nan;
xx = (xx(:,:,gap+1:end) - xx(:,:,1:end-gap)).^2;
stdMap = sqrt(nanmedian(xx,3)/0.9133);
clear xx;
stdMap(~evtSpatialMask) = nan;
stdMapGau = burst.bayesianEstimate(stdMap,T);
end