function mergingInfo = gapIdentify2D(evtLst,dF,mergingInfo,majorInfo,opts)
    [H,W,T] = size(dF);

    %% setting
%     round = 10;
%     gapSigThr = 4;
    contrastThr = opts.splitRatio;
    dw = [-1,0,1,-1,1,-1,0,1];
    dh = [-1,-1,-1,0,0,1,1,1];
    
%     %% normalize dFOrg
%     xx = (dFOrg(:,:,2:end) - dFOrg(:,:,1:end-1)).^2;
%     stdMap = sqrt(median(xx,3)/0.9133);
%     stdMapGau = double(imgaussfilt(stdMap)) + 1e-6;
%     dFOrg = dFOrg./repmat(stdMapGau,1,1,T) + opts.xBias;
    
    %%
    evtLstGrow = cell(numel(evtLst),1);
    mergingInfo.gapLst = cell(numel(evtLst),1);
    t0s = zeros(numel(evtLst),1);
    t1s = zeros(numel(evtLst),1);
    for i = 1:numel(evtLst)
        [ih0,iw0,it0] = ind2sub([H,W,T],evtLst{i});
        ihw = unique(sub2ind([H,W],ih0,iw0));
        t0s(i) = min(it0);
        t1s(i) = max(it0);
        curGrow = [];
        for k = 1:numel(dw)
           ih = max(1,min(H,ih0+dh(k)));
           iw = max(1,min(W,iw0+dw(k)));
           newAdd = sub2ind([H,W],ih,iw);
           curGrow = [curGrow;newAdd];
        end
        evtLstGrow{i} = unique(curGrow);        
    end
    
    for i = 1:numel(mergingInfo.neibLst)
        curLabel = i;
        curGrow = evtLstGrow{curLabel};
        neiLst = mergingInfo.neibLst{i};
        neiLst = neiLst(neiLst>curLabel);
        checkGap = true(numel(neiLst),1);
%         maxCurV = max(dF(evtLst{curLabel}));
        for j = 1:numel(neiLst)
            nLabel = neiLst(j);
            nGrow = evtLstGrow{nLabel};
%             maxnV = max(dF(evtLst{nLabel}));
            t0 = min(t0s(curLabel),t0s(nLabel));
            t1 = max(t1s(curLabel),t1s(nLabel));

%             ext = 5;
%             TW0 = [majorInfo{curLabel}.tPeak,majorInfo{nLabel}.tPeak];
%             t0 = max(1,min(TW0)-ext);
%             t1 = min(T,max(TW0)+ext);
            dFMean = mean(dF(:,:,t0:t1),3);
            
%             % gap
%             inter = intersect(curGrow,nGrow);
%             nBoundary = numel(inter);
%             meanIntensity = mean(dFMean(inter));
            
            % current seed
            curIhw = majorInfo{curLabel}.ihw;
            [maxV] = max(dFMean(curIhw));
            map = false(H,W);
            map(curIhw(dFMean(curIhw)>0.8*maxV)) = true;
            cc = bwconncomp(map);
            cc = cc.PixelIdxList;
            [~,reg] = max(cellfun(@numel,cc));
            tmpIhw = cc{reg};
            [curMax,curId] = max(dFMean(tmpIhw));
            curSeed = tmpIhw(curId);

            % neighbor seed
            nIhw = majorInfo{nLabel}.ihw;
            [maxV] = max(dFMean(nIhw));
            map = false(H,W);
            map(nIhw(dFMean(nIhw)>0.8*maxV)) = true;
            cc = bwconncomp(map);
            cc = cc.PixelIdxList;
            [~,reg] = max(cellfun(@numel,cc));
            tmpIhw = cc{reg};
            [nMax,nId] = max(dFMean(tmpIhw));
            nSeed = tmpIhw(nId);
            
            % new gap location
            maxSearch = min(curMax,nMax);
            minSearch = min(min(dFMean(curGrow)),min(dFMean(nGrow)));
            gap = findGap(dFMean,curSeed,nSeed,curGrow,nGrow,minSearch,maxSearch,dw,dh);
            meanIntensity = mean(dFMean(gap));
            curSeed = growMaximal(curSeed,dw,dh,H,W,curIhw,numel(gap));
            curSeedV = mean(dFMean(curSeed));
            nSeed = growMaximal(nSeed,dw,dh,H,W,nIhw,numel(gap));
            nSeedV = mean(dFMean(nSeed));
            
            
            if(contrastThr*meanIntensity<min(curSeedV,nSeedV))
                checkGap(j) = false;
            end

%             %% grow boundary
%             growRound = cell(round*2,1);
%             curGrow = inter;
%             growRound{1} = inter;
%             Range = union(evtLstGrow{curLabel},evtLstGrow{nLabel});
%             for k = 2:round*2
%                curGrow = growOneRound(curGrow,dw,dh,H,W,Range);
%                growRound{k} = curGrow;
%             end
%             
%             t0 = min(t0s(curLabel),t0s(nLabel));
%             t1 = max(t1s(curLabel),t1s(nLabel));
%             dFOrgmean = mean(dFOrg(:,:,t0:t1),3)*sqrt(t1-t0+1);
%             %% order statistic
%             for k = 1:round
%                 bg = growRound{k};
%                 fg = setdiff(growRound{2*k},bg);
%                 if(isempty(fg)||isempty(bg))
%                     continue;
%                 end
%                 fg = dFOrgmean(fg); bg = dFOrgmean(bg);
%                 L = mean(fg)-mean(bg);
%                 [mu, sigma] = ksegments_orderstatistics_fin(fg, bg);
%                 zscore = (L-mu)/sigma;
%                 if(zscore>gapSigThr)
%                     checkGap(j) = false;
%                     break;
%                 end
%             end
        end
%         mergingInfo.gapLst{i} = neiLst(~checkGap);
        neiLst = neiLst(checkGap);
        mergingInfo.neibLst{i} = neiLst;
        
    end

    if 0
       dFMean = dFMean/max(dFMean(:));
       ov = cat(3,dFMean,dFMean,dFMean);
       mask = zeros(H,W);
       mask(gap) = 1;
       ov(:,:,1) = ov(:,:,1) + mask;
       mask = zeros(H,W);
       mask(curSeed) = 1;
       ov(:,:,2) = ov(:,:,2) + mask;
       mask = zeros(H,W);
       mask(nSeed) = 1;
       ov(:,:,3) = ov(:,:,3) + mask;
       zzshow(ov)
    end
end
function pix = growOneRound(pix,dw,dh,H,W,Range)

    [ih0,iw0] = ind2sub([H,W],pix);
    curGrow = [];
    for k = 1:numel(dw)
       ih = max(1,min(H,ih0+dh(k)));
       iw = max(1,min(W,iw0+dw(k)));
       newAdd = sub2ind([H,W],ih,iw);
       curGrow = [curGrow;newAdd];
    end
    if(exist('Range','var'))
        pix = intersect(Range,curGrow);
    else
        pix = unique(curGrow);
    end
end

function pix = growMaximal(pix,dw,dh,H,W,Range,nSize)
    newAdd = pix;
%     nSize = min(nSize,10);
    while(numel(pix)<nSize)
        prePix = pix;
        [ih0,iw0] = ind2sub([H,W],newAdd);
        curGrow = [];
        for k = 1:numel(dw)
           ih = max(1,min(H,ih0+dh(k)));
           iw = max(1,min(W,iw0+dw(k)));
           newAdd = sub2ind([H,W],ih,iw);
           curGrow = [curGrow;newAdd];
        end
        curGrow = intersect(Range,curGrow);
        newAdd = setdiff(curGrow,prePix);
        if(isempty(newAdd))
            break;
        end
        pix = [pix;newAdd];
    end
end

function gap = findGap(dFMean,curSeed,nSeed,curIhw,nIhw,minSearch,maxSearch,dw,dh)
    searchThrs = minSearch:(maxSearch-minSearch)/20:maxSearch;
    % limit the gap only detected in the valid region
    valid = union(curIhw,nIhw);
    mask = false(size(dFMean));
    mask(valid) = true;
    dFMean(~mask) = minSearch-1;
    
    for i = 1:numel(searchThrs)
       BW = dFMean> searchThrs(i);
       L = bwlabel(BW);
       if(L(curSeed)~=L(nSeed))
           break;
       end
    end
    curIhw = find(L==L(curSeed));
    nIhw = find(L==L(nSeed));
    [H,W] = size(dFMean);
    gap = growUntilInterSect(curIhw,nIhw,dw,dh,H,W);
end
function inter = growUntilInterSect(curIhw,nIhw,dw,dh,H,W)

    inter = [];
    while(numel(inter)<10)
        curIhw = se.growOneRound(curIhw,dw,dh,H,W);
        nIhw = se.growOneRound(nIhw,dw,dh,H,W);
        inter = intersect(curIhw,nIhw);
    end
end