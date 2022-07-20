function [mIhw,TW,delays] = seed2Majoirty(ihw,dFVec,sz,curEvt,TW,mergedOrNot,delays)

    H = sz(1); W = sz(2); T = sz(3);
    if(~exist('mergedOrNot','var'))
        mergedOrNot = false;
    end
    if(~exist('delays','var'))
        delays = zeros(numel(ihw),1);
    end
    
    refCurve = mean(dFVec(ihw,:),1);
    TW = se.getMajorityTem(imgaussfilt(refCurve,2),TW,curEvt,[H,W,T],mergedOrNot);
    refCurve = refCurve(TW);
    L = numel(TW);
    maxtStart = T - L + 1;
    
    dw = [-1,0,1,-1,1,-1,0,1];
    dh = [-1,-1,-1,0,0,1,1,1];
    newAdd = ihw;
    delayMap = nan(H,W);
    delayMap(ihw) = delays;
    tested = false(H,W);
    tested(ihw) = true;
    
    % correlation threshold
    corSigThr = 1e-3;
    n = min(40,L);
    zscoreThr = -tinv(corSigThr,n-2);
    tmp = (zscoreThr/sqrt(n-2))^2;
    rThr = sqrt(tmp/(tmp+1));
    mIhw = ihw;    
    while(1)
        maxCor = -ones(H,W);
        delay = zeros(H,W);
        tPeaks = zeros(H,W);
        shifts = unique(delayMap(newAdd));
        candidate = [];
        for t = 1:numel(shifts)
            curShiftPixs = newAdd(delayMap(newAdd)==shifts(t));
            [ih0,iw0] = ind2sub([H,W],curShiftPixs);
            curNei = [];
            for i = 1:numel(dw) 
                ih = max(1,min(H,ih0+dh(i)));
                iw = max(1,min(W,iw0+dw(i)));
                curNei = [curNei;sub2ind([H,W],ih,iw)];
            end
            curNei = unique(curNei);
            curNei = curNei(~tested(curNei));
            candidate = [candidate;curNei];
            shift = shifts(t) + [-1:1];
            tStart = TW(1) + shift;
            tStart = tStart(tStart>0 & tStart<=maxtStart);
            for j = 1:numel(tStart)
                alignCurves = dFVec(curNei,tStart(j):tStart(j)+L-1);
                r = corr(refCurve',alignCurves')';
                select = r>maxCor(curNei);
                selectPix = curNei(select);
                maxCor(selectPix) = r(select);
                delay(selectPix) = tStart(j)-TW(1);
                [~,id] = max(imgaussfilt(dFVec(selectPix,tStart(j):tStart(j)+L-1),[1e-4,1]),[],2);
                tPeaks(selectPix) = id +  tStart(j) - 1;
            end
        end
        
        % correlation limitation
        newAdd = candidate(maxCor(candidate)>rThr);
        % peak limitation
        newAddVox = sub2ind([H*W,T],newAdd,tPeaks(newAdd));
        newAdd = newAdd(ismember(newAddVox,curEvt));
        
        % update
        mIhw = [mIhw;newAdd]; 
        tested(candidate) = true;
        delayMap(newAdd) = delay(newAdd);
        if(isempty(newAdd))
            break; 
        end    
    end
    delays = delayMap(mIhw);
end