function [spLst,cx0,tDly,neiLst,sv_spLabels,rgt00,isFail] = spgtw_MuYuProject_SVGTW_New(...
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
    cx0 = [];
    tDly = [];
    neiLst = [];
    isFail = 1;
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
stdsAll = [];

%%
neiSvLst = cell(nSv,1);
sv2DLst = cell(nSv,1);
svMap = zeros(H,W,T);
thrVec = 0.3:0.05:0.7;
medianSvDly = zeros(nSv,1);
medianDly = [];
for kk = 1:nSv
    pix = superVoxels{kk};
    [ih,iw,it] = ind2sub([H,W,T],pix);
    ihw = unique(sub2ind([H,W],ih,iw));
    svMap(superVoxels{kk}) = kk;
    sv2DLst{kk} = ihw;
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
        estPeakShift = 10;
        mtPeakTW = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T);
        mtPeakTW = intersect(mtPeakTW,mTW);
        mCurve = imgaussfilt(mean(dFVec(mIhw,:),1),2);
        dif = getMinDifInTW(mCurve,mtPeakTW,mTW);
        %
    
        curSp = mergeSmallSp(L,spSzMinLimitation,mtPeakTW,mTW,dFVec,dif);
        %% merge small super pixels
%         curSp = mergeSmallSp(L,spSzMinLimitation);
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
%     [ref,tst0,refBase,s,t,~,orgCurves0,stds] = gtw.sp2graph_alter3(dFVec,curSp,majorInfo{kk},spMap1,unique(it));
%     [ref,tst0,refBase,s,t,~,orgCurves0,stds] = gtw.sp2graph_alter4(dFVec,curSp,majorInfo{kk},spMap1,unique(it),spSz,seMap0,seSel);
    [ref,tst0,refBase,s,t,~,orgCurves0,stds,smoBases] = gtw.sp2graph_alter6(dFVec,curSp,majorInfo{kk},spMap1,unique(it),spSz,seMap0,seSel,smoBase);
    nSp = numel(curSp);
    sz = cellfun(@numel,curSp);
    s2 = stds;
    
    %% GTW
    [ ss,ee,gInfo ] = gtw.buildGTWGraph( ref, tst0, s, t, smoBases, maxStp, s2);
    [~, labels1] = aoIBFS.graphCutMex(ss,ee);
    path0 = gtw.label2path4Aosokin( labels1, ee, ss, gInfo );
    
    %% warp curves
    cx0 = gtw.warpRef2Tst_alter(path0,refBase/max(refBase(:)),[H,W,numel(refBase)]);
    cx = [cx;cx0];
    tst = [tst;tst0];
    stdsAll = [stdsAll;stds];
    orgCurves = [orgCurves;orgCurves0];
    %% time to achieve different levels for each seed
    
    tAch = nan(nSp,numel(thrVec));
    for nn=1:nSp
        x = cx0(nn,:);
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
    tDly0 = mean(tAch,2);
    tDly(spLabels) = tDly0;
    medianSvDly(kk) = median(tDly0);
%     medianDly(spLabels) = medianSvDly(kk);
end
neiSvLst = getNeiLst(svMap,sv_spLabels,sv2DLst,size(dF));
smoMedianSvDly = zeros(nSv,1);
for k = 1:nSv
    neiSv = neiSvLst{k};
    neiSv = [neiSvLst{k},k];
    smoMedianSvDly(k) = round(mean(medianSvDly(neiSv)));
    medianDly(sv_spLabels{k}) = smoMedianSvDly(k);
end

% tDly = tDly';
thrVec = 0.5:0.05:0.95;
m0Msk = seMap0==seSel;
validMap = sum(m0Msk,3)>0;
dFip = dF;
dFip(seMap0~=seSel) = nan;
dFip = gtw.imputeMov_Fast(dFip,validMap);
dFAvg = max(dFip,[],3);
[~,seedIn] = max(dFAvg);
[ih0,iw0] = find(validMap>0);
gapSeed = max(ceil(max(max(ih0)-min(ih0),max(iw0)-min(iw0))/opts.gtwGapSeedRatio),opts.gtwGapSeedMin);
[refBase,baseDly] = getRefBase(dFip,validMap,seedIn,gapSeed,thrVec);
[refAll,tstShift,refBase] = shiftTstCurve(tst,medianDly,refBase,baseDly);
[neiLst,sAll,tAll] = getNeiLst(spMap,sv_spLabels,spLst,size(dF));
smoBases = ones(nSpAll,1)*smoBase;
neiNum = cellfun(@numel,neiLst);
smoBases(neiNum<3) = smoBase*2;

%% GTW
[ ss,ee,gInfo ] = gtw.buildGTWGraph( refAll, tstShift, sAll, tAll, smoBases, maxStp, stdsAll);
[~, labels1] = aoIBFS.graphCutMex(ss,ee);
path0 = gtw.label2path4Aosokin( labels1, ee, ss, gInfo );
cx0 = gtw.warpRef2Tst_alter(path0,refBase,[H,W,numel(refBase)]);
tAch = nan(nSpAll,numel(thrVec));
for nn=1:nSpAll
    x = cx0(nn,:);
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
tDly = mean(tAch,2) + medianDly' - baseDly;

% %% make up false sp
% for ii=1:numel(falseSp)
%     curLabel = falseSp(ii);
%     neib0 = neiLst{curLabel};
%     neib0 = setdiff(neib0,falseSp);
%     tDly(curLabel) = mean(tDly(neib0));
% end
 % [svLabel] = burst.riseMap2evtMuYu_Alter2(spLst,tDly,neiLst,sv_spLabels,6,6,[H,W]);
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
function [refBase,medianDly] = getRefBase(dFip,validMap,seedIn,gapSeedHW,thrVec)
    [H,W,T] = size(dFip);
    [ih,iw] = ind2sub([H,W],seedIn);
    rgh = max(ih-gapSeedHW,1):min(ih+gapSeedHW,H);
    rgw = max(iw-gapSeedHW,1):min(iw+gapSeedHW,W);
    df00Vec = reshape(dFip(rgh,rgw,:),[],T);
    vm00 = validMap(rgh,rgw);
    df00Vec = df00Vec(vm00>0,:);
    refBase = nanmean(df00Vec,1);
    
    [~,tPeak] = max(imgaussfilt(refBase,2));
    ext = 5; ts = 1;
    validPart = false(T,1);
    validPart(ts:min(tPeak+ext,T)) = true;
    %% make it non-decreasing
    smo = 1;
    while(1)
        curveSmo = imgaussfilt(refBase,smo);
        left = false(T,1);
        left(1) = true;
        left(2:end) = curveSmo(2:end)>=curveSmo(1:end-1);
        right = false(T,1);
        right(end) = true;
        right(1:end-1) = curveSmo(1:end-1)>=curveSmo(2:end);
        maxima = left & right & validPart;
        maxima = find(maxima);
        if(numel(maxima)<=1 || smo>=5)
           break;
        end
        smo = smo + 1;
    end
    refBase = curveSmo;
    %%
    peakTW = ts:min(tPeak + ext,T);
    [~,tPeak] = max(curveSmo(peakTW));
    tPeak = ts + tPeak - 1;
    refBase(1:ts) = refBase(ts);
    refBase(tPeak:end) = refBase(tPeak);
    
    refBase = refBase - refBase(1);
    refBase = refBase/refBase(tPeak);

    %% non decreasing
    preV = 1;
    preT = tPeak;
    for t = tPeak:-1:ts
       if(refBase(t)<preV) 
           refBase(t:preT) = refBase(t) + (preV-refBase(t))/(preT-t)*[0:preT-t];
           preV = refBase(t);
           preT = t;
       end
    end
    
    tAchRef = nan(1,numel(thrVec));
    x = refBase;
    [~,t0] = max(x);
    x = x(1:t0);
    for ii=1:numel(thrVec)
        t1 = find(x>=thrVec(ii),1);
        if isempty(t1)
            t1 = t0;
        end
        tAchRef(ii) = t1;
    end
    medianDly = round(mean(tAchRef));
end
function [ref,tstShift,refBase] = shiftTstCurve(tst,medianDly,refBase,baseDly)
    nSp = size(tst,1);
    extPre = max(max(medianDly-baseDly),0);
    extPost = max(max(baseDly-medianDly),0);
    refBase = [zeros(1,extPre),refBase,ones(1,extPost)];
    tst = [zeros(nSp,extPre),tst,ones(nSp,extPost).*max(tst,[],2)];
    tstShift = zeros(size(tst));
    ref = zeros(size(tst));
    T = size(tstShift,2);
    tss = T;
    tee = 0;
    ext = 10;
    for i = 1:nSp
        shift = medianDly(i) - baseDly;
        vMax = max(tst(i,:));
        if(shift>=0)
            tstShift(i,1:end-shift) = tst(i,shift+1:end);
            tstShift(i,end-shift+1:end) = vMax*ones(1,shift);
        else
            shift = -shift;
            tstShift(i,1:shift) = zeros(1,shift);
            tstShift(i,shift+1:end) = tst(i,1:end-shift);
        end
        ref(i,:) = refBase*vMax;
        t00 = find(tstShift(i,:)~=0,1)-1;
        t11 = find(tstShift(i,:)<vMax,1,'last'); 
        if(isempty(t11)|| isempty(t11)) 
            continue; 
        end
        tss = min(t00,tss);
        tee = max(t11,tee);
    end
    tss = max(tss - ext,1);
    tee = min(tee + ext,T);
    
    tstShift = tstShift(:,tss:tee);
    ref = ref(:,tss:tee);
    refBase = refBase(tss:tee);
end
function [neiLst,sAll,tAll] = getNeiLst(spMap,sv_spLabels,spLst,sz)
    H = sz(1); W = sz(2); T = sz(3); nSpAll = numel(spLst);
    
    dh = [-1 0 1 -1 0 1 -1 0 1];
    dw = [-1 -1 -1 0 0 0 1 1 1];
    neiLst = cell(nSpAll,1);
    spMap = reshape(spMap,[],T);
    sAll = [];
    tAll = [];

    % sp_svLabel = zeros(numel(nSpAll),1);
    % for i = 1:numel(sv_spLabels)
    %     curLabels = sv_spLabels{i};
    %     sp_svLabel(curLabels) = i;
    % end
    % 
    % for ii = 1:nSpAll
    %     svLabel = sp_svLabel(ii);
    %     sameSv = sv_spLabels{svLabel};
    %     sp0 = spLst{ii};
    %     [ih0,iw0] = ind2sub([H,W],sp0);
    %     neib0 = [];
    %     for jj=1:numel(dh)  % find neighbors in eight directions
    %         ih = max(min(H,ih0+dh(jj)),1);
    %         iw = max(min(W,iw0+dw(jj)),1);
    %         pixShift = sub2ind([H,W],ih,iw);
    %         neib0 = union(neib0,spMap(pixShift,:));
    %     end
    %     neib0 = neib0(neib0>ii);
    %     neib0 = intersect(neib0,sameSv);
    %     for j = 1:numel(neib0)
    %         s0 = ii;
    %         t0 = neib0(j);
    %         neiLst{s0} = union(neiLst{s0},t0);
    %         neiLst{t0} = union(neiLst{t0},s0);
    %     end
    % end
    % spMap = reshape(spMap,[H,W,T]);
    % 
    % for t = 1:T
    %     tmp = spMap(:,:,t);
    %     mask = tmp>0;
    %     mask = imerode(mask,strel('square',3));
    %     tmp(~mask) = 0;
    %     spMap(:,:,t) = tmp;
    % end
    % spMap = reshape(spMap,[],T);


    % formal search neighbors
    for ii=1:nSpAll
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
        for j = 1:numel(neib0)
            s0 = ii;
            t0 = neib0(j);
            neiLst{s0} = union(neiLst{s0},t0);
            neiLst{t0} = union(neiLst{t0},s0);
        end
        sAll = [sAll;ones(numel(neib0),1)*ii];
        tAll = [tAll;neib0];
    end
end