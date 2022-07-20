function [spLst,cx,tDly,neiLst,sv_spLabels,rgt00,isFail] = spgtw_MuYuProject_DouDi(...
    dF,seMap,seSel,smoBase,maxStp,cDelay,spSz,spT,superVoxels,majorInfo,opts)
% spgtw super pixel GTW 
% make one burst to super pixels and run gtw
% tic;
if ~isfield(opts,'gtwGapSeedMin') || ~isfield(opts,'gtwGapSeedRatio')
    opts.gtwGapSeedRatio = 4;
    opts.gtwGapSeedMin = 5;
end

%% setting
[H,W,T] = size(dF);
isFail = 0;
maxStp = max(min(maxStp,round(T/2)),1);

%% super pixels
m0Msk = seMap==seSel;
validMap = sum(m0Msk,3)>0;
dFip = dF;
dFip(seMap~=seSel) = nan;
dFip = gtw.imputeMov_Fast(dFip,validMap);
dFAvg = max(dFip,[],3);
%% super pixel
spSzMinLimitation = spSz/4;
nSpEstimate = max(1,round(H*W/spSz));

ts = T;
te = 1;
for kk = 1:numel(superVoxels)
    [~,~,it] = ind2sub([H,W,T],superVoxels{kk});
    ts = min(ts,min(it));
    te = max(te,max(it));
end
dFAvg2 = max(dF(:,:,ts:te),[],3);
[L] = superpixels(dFAvg2,nSpEstimate,'Compactness',20);
% L(~validMap) = 0;
spLst = label2idx(L);

%% get super pixel
if sum(validMap(:))< spSzMinLimitation || numel(spLst)<2
    spLst = {find(validMap>0)};
    cx = [];
    tDly = [];
    neiLst = [];
    rgt00 = 1:T;
    isFail = 1;
    sv_spLabels = [];
    return
end

% %% connectivity
% nSp = numel(spLst);
% nAdd = nSp + 1;
% for i = 1:nSp
%     pix = spLst{i};
%    if(~isempty(pix)) 
%        cMap = false(H,W);
%        cMap(pix) = true;
%        cc = bwconncomp(cMap);
%        if(cc.NumObjects>1)
%            comps = cc.PixelIdxList;
%            for j = 1:numel(comps)
%                if(j==1)
%                    spLst{i} = comps{j};
%                else
%                    spLst{nAdd} = comps{j};
%                    L(comps{j}) = nAdd;
%                    nAdd = nAdd+1;
%                end
%            end
%        end
%    end
% end

%% 
dFVec = reshape(dF,[],T);
spLst = cell(0);
nSp = 0;
sv_spLabels = cell(numel(superVoxels),1);
for k = 1:numel(superVoxels)
    pix = superVoxels{k};
    [ih,iw,~] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    L0 = zeros(H,W);
    L0(ihw) = L(ihw);
    % add check signal
    
    %% connectivity
%     [L0] = spConn(L0);
        
    mIhw = majorInfo{k}.ihw;
    mTW = majorInfo{k}.TW;
    mtPeak = majorInfo{k}.tPeak;
    estPeakShift = 10;
    mtPeakTW = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T);
    mtPeakTW = intersect(mtPeakTW,mTW);
    mCurve = imgaussfilt(mean(dFVec(mIhw,:),1),2);
    dif = getMinDifInTW(mCurve,mtPeakTW,mTW);
    %
    
    curSp = mergeSmallSp(L0,spSzMinLimitation,mtPeakTW,mTW,dFVec,dif);
    spLabels = nSp + [1:numel(curSp)];
    nSp = nSp + numel(curSp);
    spLst(spLabels) = curSp;
    sv_spLabels{k} = spLabels;
end

%%
[~,seedIn] = max(dFAvg);

% signal part
idx0 = find(m0Msk>0);
[~,~,it0] = ind2sub(size(m0Msk),idx0);
rgt00 = max(min(it0)-5,1):min(max(it0)+5,T);
dat = double(dFip(:,:,rgt00));

fprintf('Node %d\n',numel(spLst));

%% alignment
% graph
[ih0,iw0] = find(validMap>0);
gapSeed = max(ceil(max(max(ih0)-min(ih0),max(iw0)-min(iw0))/opts.gtwGapSeedRatio),opts.gtwGapSeedMin);
% [ref,tst,refBase,s,t,idxGood,orgCurves] = gtw.sp2graph_alter2(dat,dF(:,:,rgt00),m0Msk(:,:,rgt00),validMap,spLst,seedIn,gapSeed);
[ref,tst,refBase,s,t,idxGood,orgCurves,s2] = gtw.sp2graph_alterDouDi(dat,dF,m0Msk,rgt00,validMap,spLst,seedIn,gapSeed,majorInfo,sv_spLabels);
% toc;
% gtw
% tic;
spLst = spLst(idxGood);
nSp = numel(spLst);
sz = cellfun(@numel,spLst);
s2 = 1./sqrt(sz);
if numel(spLst)>3 && numel(refBase)>5
%     tic
    [ ss,ee,gInfo ] = gtw.buildGTWGraph( ref, tst, s, t, smoBase, maxStp, s2);
    [~, labels1] = aoIBFS.graphCutMex(ss,ee);
    path0 = gtw.label2path4Aosokin( labels1, ee, ss, gInfo );
%     t00 = toc;
%     if numel(spLst)>1000
%         fprintf('Time %fs\n',t00);
%     end
else
    [nPix,nTps] = size(tst);
    path0 = cell(nPix,1);
    rg = (0:nTps)';
    p0 = [rg,rg,rg+1,rg+1];
    for ii=1:nPix
        path0{ii} = p0;
    end
end

% warped curves
cx = gtw.warpRef2Tst_alter(path0,refBase/max(refBase(:)),[H,W,numel(refBase)]);
% toc;
% time to achieve different levels for each seed
nSp = numel(spLst);
thrVec = 0.5:0.05:0.95;
tAch = nan(nSp,numel(thrVec));
for nn=1:nSp
    x = cx(nn,:);
    [~,t0] = max(x);
    x = x(1:t0);
    for ii=1:numel(thrVec)
        t1 = find(x>=thrVec(ii),1);
        if isempty(t1)
            t1 = t0;
        end
        tAch(nn,ii) = t1;
    end
end
tDly = mean(tAch,2);

% obtain neighbor list
neiLst = cell(nSp,1);
nPair = numel(s);
for nn=1:nPair
    s0 = s(nn);
    t0 = t(nn);
    neiLst{s0} = [neiLst{s0};t0];
    neiLst{t0} = [neiLst{t0};s0];
    
%     d0 = tAch(s0,:)-tAch(t0,:);  % negative is earlier
%     d0 = sum(d0)/numel(thrVec);
%     distMat(s0,t0) = d0;
%     distMat(t0,s0) = -d0;
end

% [svLabel] = burst.riseMap2evtMuYu_Alter2(spLst,tDly,neiLst,sv_spLabels,6,6,[H,W]);
end

function spLst = mergeSmallSp(L,spSzMinLimitation,mtPeakTW,mTW,dFVec,dif)
    spLst = label2idx(L);
    [H,W] = size(L);
    dh = [-1,-1,-1,0,0,1,1,1];
    dw = [-1,0,1,-1,1,-1,0,1];
    
    isolatedRegion = [];
    sz = cellfun(@numel,spLst);
    checkRegions = find(sz>0);
    % merge small pixels
    while(1)
        changeHappen = false;
        for i = 1:numel(checkRegions)
            curLabel = checkRegions(i);
            pix = spLst{curLabel};
            curCurve = mean(dFVec(pix,:),1);
            s0 = mean((curCurve(2:end)-curCurve(1:end-1)).^2)/2;
            curDif = getMinDifInTW(imgaussfilt(curCurve,2),mtPeakTW,mTW);
            noSignal = curDif<max(dif/3,5*s0);
            curSz = numel(spLst{curLabel});
            if(curSz<spSzMinLimitation || noSignal)
                [ih,iw] = ind2sub([H,W],pix);
               for k = 1:numel(dw)
                   ih1 = min(H,max(1,ih+dh(k)));
                   iw1 = min(W,max(1,iw+dw(k)));
                   pixShift = sub2ind([H,W],ih1,iw1);
                   neiL = setdiff(L(pixShift),[0,curLabel]);
                   if(~isempty(neiL))
                      break; 
                   end
               end
               if(isempty(neiL))
                  isolatedRegion = [isolatedRegion,curLabel];
                  continue;
               end
               neiL = neiL(1);
               % update
               spLst{neiL} = [spLst{neiL};pix];
               L(pix) = neiL;
               spLst{curLabel} = [];
               changeHappen = true;
            end
        end
        checkRegions = setdiff(checkRegions,isolatedRegion);
        if(~changeHappen)
           break; 
        end
    end 
    sz = cellfun(@numel,spLst);
    spLst = spLst(sz>0);
end
function dif = getMinDifInTW(mCurve,mtPeakTW,mTW)
    [vMax,tPeak] = max(mCurve(mtPeakTW));
    tPeak = tPeak + min(mtPeakTW) - 1;
    t0 = min(mTW); t1 = max(mTW);
    vMin1 = min(mCurve(t0:tPeak));
    vMin2 = min(mCurve(tPeak:t1));
    dif = min(vMax-vMin1,vMax-vMin2);
end
function [L] = spConn(L)
	[H,W] = size(L);
    spLst00 = label2idx(L);
    nSp = numel(spLst00);
    nAdd = nSp + 1;
    for i = 1:nSp
        pix = spLst00{i};
       if(~isempty(pix)) 
           cMap = false(H,W);
           cMap(pix) = true;
           cc = bwconncomp(cMap);
           if(cc.NumObjects>1)
               comps = cc.PixelIdxList;
               for j = 1:numel(comps)
                   if(j==1)
                       spLst00{i} = comps{j};
                   else
                       spLst00{nAdd} = comps{j};
                       L(comps{j}) = nAdd;
                       nAdd = nAdd+1;
                   end
               end
           end
       end
    end
end