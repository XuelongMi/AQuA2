function [spLst,cx,tDly,neiLst,sv_spLabels,rgt00,isFail,riseDur] = spgtw_MuYuProject_GlobalSVGTW(...
    dF,seMap0,seSel,smoBase,maxStp,cDelay,spSz,spT,superVoxels,majorInfo,opts)
% spgtw super pixel GTW 
% make one burst to super pixels and run gtw

if ~isfield(opts,'gtwGapSeedMin') || ~isfield(opts,'gtwGapSeedRatio')
    opts.gtwGapSeedRatio = 4;
    opts.gtwGapSeedMin = 5;
end
% spSz = 500;
isFail = 0;
[H,W,T] = size(dF);
nSv = numel(superVoxels);
spSzMinLimitation = spSz/4;
nSpEstimate = max(1,round(H*W/spSz));
rgt00 = 1:T;

if nSv == 1
    pix = superVoxels{1};
    [ih,iw,it] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    spLst = {ihw};
    cx = [];
    tDly = [];
    neiLst = [];
    isFail = 1;
    riseDur = 1;
    sv_spLabels = [];
    return
end

ts = T;
te = 1;
for kk = 1:nSv
    [~,~,it] = ind2sub([H,W,T],superVoxels{kk});
    ts = min(ts,min(it));
    te = max(te,max(it));
end
dFAvg = max(dF(:,:,ts:te),[],3);
L_All = superpixels(dFAvg,nSpEstimate,'Compactness',20);

nSpAll = 0;
spLst = cell(0);
sv_spLabels = cell(nSv,1);
tDly = [];
dFVec = reshape(dF,[],T);
spMap = zeros(H,W,T);
falseSp = [];
cx = [];
orgCurves = [];
tst = [];
riseDur = [];
stds = [];

for kk = 1:nSv
    pix = superVoxels{kk};
    [ih,iw,it] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    maxStp = max(min(maxStp,round(T/2)),1);
    
    if(isempty(majorInfo{kk}.ihw))
        curSp = cell(1,1);
        curSp{1} = ihw;
        falseSp = [falseSp;1+nSpAll];
    else
        L = zeros(H,W);
        L(ihw) = L_All(ihw);
        %% connectivity
        [L] = spConn(L);
        % add check signal
        mIhw = majorInfo{kk}.ihw;
        mTW = majorInfo{kk}.TW;
        mtPeak = majorInfo{kk}.tPeak;
        estPeakShift = max(5,round(numel(mTW)/6));
        mtPeakTW = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T);
        mtPeakTW = intersect(mtPeakTW,mTW);
        mCurve = imgaussfilt(mean(dFVec(mIhw,:),1),2);
        dif = getMinDifInTW(mCurve,mtPeakTW,mTW);
        curSp = mergeSmallSp(L,spSzMinLimitation,mtPeakTW,mTW,dFVec,dif);
    end
    
    %% update spMap
    spMap0 = zeros(H*W,T);
    for i = 1:numel(curSp)
        spMap0(curSp{i},:) = i;
    end
    spMap0 = reshape(spMap0,[H,W,T]);
    spMap1 = zeros(H,W,T);
    spMap1(pix) = spMap0(pix);
    clear spMap0;
    spMap(pix) = spMap1(pix) + nSpAll;
    spLabels = nSpAll + [1:numel(curSp)];
    nSpAll = nSpAll + numel(curSp);
    spLst(spLabels) = curSp;
    sv_spLabels{kk} = spLabels;
    
    %% extract signals - TODO modify. major TW should be changed.
    [ref,tst0,refBase,s,t,~,orgCurves0,std0,smoBases,riseDur0] = gtw.sp2graph_alter6(dFVec,curSp,majorInfo{kk},spMap1,unique(it),spSz,seMap0,seSel,smoBase);
    nSp = numel(curSp);
    stds = [stds;std0];
    tst = [tst;tst0];
    riseDur = [riseDur;riseDur0];
    orgCurves = [orgCurves;orgCurves0];
end
nSp = size(tst,1);
normalTst = tst./repmat(max(tst,[],2),1,T);
mintst = min(normalTst,[],1);
maxtst = max(normalTst,[],1);
[refBase,maxShift] = getRefMiddle(mintst,maxtst);
% cut non-necessary part
ext = 5;
t00 = find(maxtst>0,1); t00 = max(1,t00-ext);
t11 = find(mintst<1,1,'last'); t11 = min(T,t11 + ext);
refBase = refBase(t00:t11);
tst = tst(:,t00:t11);
ref = repmat(refBase,nSp,1).*repmat(max(tst,[],2),1,numel(refBase));

% formal search neighbors
dh = [-1 0 1 -1 0 1 -1 0 1];
dw = [-1 -1 -1 0 0 0 1 1 1];
neiLst = cell(nSpAll,1);
spMap = reshape(spMap,[],T);
nPair = 0;
s = nan(nSp,1);
t = nan(nSp,1);
for ii=1:nSp
    sp0 = spLst{ii};
    [ih0,iw0] = ind2sub([H,W],sp0);
    neib0 = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H,ih0+dh(jj)),1);
        iw = max(min(W,iw0+dw(jj)),1);
        pixShift = sub2ind([H,W],ih,iw);
        neib0 = union(neib0,spMap(pixShift,:));
    end
    neib0 = neib0(neib0>ii);  
    ids = nPair+[1:numel(neib0)];
    nPair = nPair + numel(neib0);
    s(ids) = ii;
    t(ids) = neib0;
    for j = 1:numel(neib0)
        s0 = ii;
        t0 = neib0(j);
        neiLst{s0} = union(neiLst{s0},t0);
        neiLst{t0} = union(neiLst{t0},s0);
    end
end
s = s(~isnan(s));
t = t(~isnan(t));
smoBases = ones(nSp,1)*smoBase;
neiSz = cellfun(@numel,neiLst);
smoBases(neiSz<3) = 2*smoBase;

%% GTW
% maxShift = max(maxShift,opts.maxStp);
[ ss,ee,gInfo ] = gtw.buildGTWGraph( ref, tst, s, t, smoBases, maxShift, stds);
[~, labels1] = aoIBFS.graphCutMex(ss,ee);
path0 = gtw.label2path4Aosokin( labels1, ee, ss, gInfo );

%% warp curves
cx = gtw.warpRef2Tst_alter(path0,refBase/max(refBase(:)),[H,W,numel(refBase)]);

%% time to achieve different levels for each seed
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
tDly = tDly + t00 - 1;

% %% make up false sp
% for ii=1:numel(falseSp)
%     curLabel = falseSp(ii);
%     neib0 = neiLst{curLabel};
%     neib0 = setdiff(neib0,falseSp);
%     tDly(curLabel) = mean(tDly(neib0));
% end
 % [svLabel] = burst.riseMap2evtMuYu_Alter3(spLst,tDly,neiLst,sv_spLabels,6,4,[H,W],spSz,riseDur);
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
            noSignal = curDif<5*s0;
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
function [refCurve,maxShift] = getRefMiddle(mintst,maxtst)
    thrVec = 0:0.01:0.99;
    refCurve = zeros(size(mintst));
    maxShift = 0;
    for ii=1:numel(thrVec)
        t1 = find(mintst>thrVec(ii),1);
        t2 = find(maxtst>thrVec(ii),1);
        midT = floor((t1+t2)/2);
        maxShift = max([maxShift,abs(t1-midT),abs(t2-midT)]);
        refCurve(midT) = thrVec(ii)+1e-4;
    end
    % start
    t1 = find(mintst>0,1);
    t2 = find(maxtst>0,1);
    midT = floor((t1+t2)/2) -1;
    maxShift = max([maxShift,abs(t1-midT),abs(t2-midT)]);
    refCurve(1:midT) = 0;
    % end
    t1 = find(mintst==1,1);
    t2 = find(maxtst==1,1);
    midT = floor((t1+t2)/2);
    maxShift = max([maxShift,abs(t1-midT),abs(t2-midT)]);
    refCurve(midT:end) = 1;
    preV = 0;
    preT = 1;
    for t = 2:numel(mintst)
       if(refCurve(t)>=preV) 
           refCurve(preT:t) = refCurve(preT) + (refCurve(t)-preV)/(t-preT)*[0:t-preT];
           preV = refCurve(t);
           preT = t;
       end
    end
    % tolerance
    maxShift = maxShift + 2;
end