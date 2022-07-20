function [outputArg1,outputArg2] = curSp(lblMap,dF)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[H,W,T] = size(dF);
dFVec = reshape(dF,[],T); clear dat;
spLst = label2idx(lblMapS);
nSp = numel(spLst);
for i = 1:numel(spLst)
    [ih,iw,it] = ind2sub([H,W,T],spLst{i});
    ihw = unique(sub2ind([H,W],ih,iw));
    t0 = min(it);
    t1 = max(it);
    curve = mean(dFVec(ihw,t0:t1),1);
    
    [peakValue,tPeak] = max(curve);
    tw = cell(0);
    
    while()
    
    
    
    
    
end




end

