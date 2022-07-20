function [spLst,tDly,dlyMap,majorityMap,scl] = spgtw3(...
    dF,seMap0,seSel,smoBase,spSz,superVoxels,majorInfo)
% spgtw super pixel GTW 
% make one burst to super pixels and run gtw;
[H,W,T] = size(dF);
nSv = numel(superVoxels);
spSzMinLimitation = spSz/4;
nSpEstimate = max(1,round(H*W/spSz));

ts = T;
te = 1;
for kk = 1:nSv
    [~,~,it] = ind2sub([H,W,T],superVoxels{kk});
    ts = min(ts,min(it));
    te = max(te,max(it));
end
dFAvg = max(dF(:,:,ts:te),[],3);
spMap = superpixels(dFAvg,nSpEstimate,'Compactness',20);
% dF(seMap0~=seSel) = 3;
% dF = dF-3;
% dF = gtw.imputeMov_Fast(dF,sum(seMap0==seSel,3)>0);
dFVec = reshape(dF,[],T);
% clear dF;

%% super pixels
% draw the majority map of this super event. 
% majority may have overlap, thus in temporal order. Earlier is better
majorityMap = nan(H,W);
majorityMap(sum(seMap0==seSel,3)>0) = 0;
majorityMap0 = zeros(H,W);
tRising = zeros(nSv,1);     
for kk = 1:nSv
    tRising(kk) = majorInfo{kk}.tPeak;
end
[~,ids] = sort(tRising,'descend');
for kk = 1:nSv
    curId = ids(kk);
    majorityMap(majorInfo{curId}.ihw) = curId;
    [ih,iw,~] = ind2sub([H,W,T],superVoxels{curId});
    majorityMap0(unique(sub2ind([H,W],ih,iw))) = kk;
end

% divide regions into super voxels, let one super pixel could only have one
% label
spLst = label2idx(spMap);
nSp = numel(spLst);
for k = 1:nSp
    pix = spLst{k};
    labels = unique(majorityMap(pix));
    labels = labels(~isnan(labels));
    spMap(pix) = 0;
    if (~isempty(labels))
        for i = 1:numel(labels)
            curLabel = labels(i);
            pix0 = pix(majorityMap(pix)==curLabel);
            if(i==1)
                spMap(pix0) = k;
            else
                nSp = nSp + 1;
                spMap(pix0) = nSp;
            end
        end
    end
end
spLst = label2idx(spMap);
sz = cellfun(@numel,spLst);
spLst = spLst(sz>0);
nSp = numel(spLst);
% update
for k = 1:nSp
    spMap(spLst{k}) = k;
end
% refine: merge small super pixel
isolated = false(1,nSp);
dh = [-1,-1,-1,0,0,1,1,1];
dw = [-1,0,1,-1,1,-1,0,1];
while(1)
    sz = cellfun(@numel,spLst);
    needMerge = find(sz<spSzMinLimitation & sz > 0 & ~isolated);
    if (isempty(needMerge))
        break;
    end
    for i = 1:numel(needMerge)
        curLabel = needMerge(i);
        pix = spLst{curLabel};
        [ih,iw] = ind2sub([H,W],pix);
        pixShift = [];
        for k = 1:numel(dw)
           ih1 = min(H,max(1,ih+dh(k)));
           iw1 = min(W,max(1,iw+dw(k)));
           pixShift = [pixShift;sub2ind([H,W],ih1,iw1)];
        end
        pixShift = setdiff(pixShift,pix);
        pixShift = pixShift(majorityMap(pixShift)==majorityMap(pix(1)));
        if(isempty(pixShift))
            isolated(curLabel) = true;
        else
           % merge
           mergeLabel = spMap(pixShift(1));
           spLst{mergeLabel} = [spLst{mergeLabel};pix];
           spLst{curLabel} = [];
        end

    end
end
sz = cellfun(@numel,spLst);
spLst = spLst(sz>0);
spMap = zeros(H,W);
for i = 1:numel(spLst)
    spMap(spLst{i}) = i;
end
nSp = numel(spLst);


%% extract signals
% temporal downsample, maximum 30 time points
scl = max(1,floor(T/30));
nonvalidMap = reshape(seMap0>0 & seMap0~=seSel,[],T);
if (scl>1)
    T0 = ceil(T/scl);
    dFDS = zeros(H*W,T0);
    nonvalidMapDS = false(H*W,T0);
    for t = 1:T0
       t0 = (t-1)*scl + 1;
       t1 = min(T,t*scl);
       dFDS(:,t) = mean(dFVec(:,t0:t1),2);
       nonvalidMapDS(:,t) = max(nonvalidMap(:,t0:t1),[],2);
    end
    dFVec = dFDS;
    nonvalidMap = nonvalidMapDS;
    clear dFDS validMapDS;
    T = T0;
end
% pick the largest majority's curve as reference curve
refBase = 0:1/(T-1):1;
% refIhw = [];
% for k = 1:nSv
%     ihw = majorInfo{k}.ihw;    
%     if(numel(ihw)>numel(refIhw))
%        refIhw = ihw;
%     end
% end
% refBase = nanmean(dFVec(refIhw,:),1);
% % refine refBase
% refBase = imgaussfilt(refBase,2);
% TW = max(~nonvalidMap(refIhw,:),[],1);
% [~,tPeak] = max(refBase(TW));
% tPeak = tPeak + find(TW,1) - 1;
% for t = tPeak-1:-1:1
%    refBase(t) = min(refBase(t+1),refBase(t));
% end
% for t = tPeak+1:T
%    refBase(t) = min(refBase(t-1),refBase(t));
% end
% refBase = refBase - min(refBase);
% refBase = refBase/max(refBase);
% refBase(tPeak:end) = 1;
% aligned curves
ref = zeros(nSp,T);
tst = zeros(nSp,T);
ext = 5;
for k = 1:nSp
    ihw = spLst{k};
    id = mode(majorityMap0(ihw));
    TW = unique(ceil((majorInfo{id}.TW)/scl));
%     prePeak = ceil((majorInfo{id}.tPeak + 1)/scl);
%     pTW = max(prePeak-ext,1):min(prePeak + ext,T);
    curve = mean(dFVec(ihw,:),1)*sqrt(numel(ihw));
    curve = imgaussfilt(curve,2);
    [maxV,tPeak] = max(curve(TW));
    tPeak = tPeak + TW(1) - 1;
    ts1 = max(min(TW(1)-ext),1);
    [minV,ts] = min(curve(ts1:tPeak));
    ts = ts1 + ts - 1;
    curve(1:ts) = minV;
    curve(tPeak:end) = maxV;
    curve = curve - minV;
    curve = curve - min(curve);
    tst(k,:) = curve;
    ref(k,:) = refBase * (maxV-minV);
end
% clear dFVec; clear seMap0;

% get neighbor relation
dh = [-1 0 1 -1 0 1 -1 0 1];
dw = [-1 -1 -1 0 0 0 1 1 1];
Gij = zeros(nSp*10,2);
nPair = 0;
for i = 1:nSp
    ihw = spLst{i};
    [ih0,iw0] = ind2sub([H,W],ihw);
    neib0 = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H,ih0+dh(jj)),1);
        iw = max(min(W,iw0+dw(jj)),1);
        pixShift = sub2ind([H,W],ih,iw);
        neib0 = union(neib0,spMap(pixShift));
    end
    neib0 = neib0(neib0>i);  
    ids = nPair+[1:numel(neib0)];
    nPair = nPair + numel(neib0);
    Gij(ids,1) = i;
    Gij(ids,2) = neib0;
end
Gij = Gij(1:nPair,:);

%% GTW: BILCO calculate matching
distMatrix = zeros(T,T,size(tst,1));
for i = 1:nSp
   distMatrix(:,:,i) = (ref(i,:)-tst(i,:)').^2;
end
if(smoBase==0 || isempty(Gij))
    midPoints = zeros(nSp,T-1);
    for i = 1:nSp
        midPoints(i,:) = DTW_Edge_input(distMatrix(:,:,i));
    end
else
    initialCut0 = DTW_Edge_input(mean(distMatrix,3))+1;
    initialCut = repmat(initialCut0,size(tst,1),1);
    clear distMatrix;
    midPoints = BILCO(ref,tst,Gij,smoBase,initialCut);
end
paths = cell(size(tst,1),1);
for i = 1:size(midPoints,1)
   paths{i} = midPoint2path(midPoints(i,:),T,T);
end
%% warping curves
cx = gtw.warpRef2Tst_alter(paths,refBase,[H,W,T]);

%% delay time
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

dlyMap = nan(H,W);
for i = 1:nSp
    dlyMap(spLst{i}) = tDly(i);
end
end