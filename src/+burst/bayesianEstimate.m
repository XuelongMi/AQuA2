function [stdMapGau] = bayesianEstimate(stdMap,T)
    mu0 = nanmean(stdMap(:));
    sig0square = nanstd(stdMap(:)).^2;
    stdMap(isnan(stdMap)) = mu0;
    sig1square = stdMap.^2/2/T;
    stdMap_afterDeal = (stdMap*sig0square + mu0*sig1square)./(sig0square+sig1square);
    stdMapGau = double(imgaussfilt(stdMap_afterDeal)) + 1e-6;
    stdMapGau(stdMapGau<max(stdMapGau(:))/30) = max(stdMapGau(:))/30;
end

