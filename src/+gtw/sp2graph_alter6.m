function [ref,tst,refBase,s,t,idxGood,orgCurves,stds,smoBases,riseDur] = sp2graph_alter6(dFVec,spLst,major,spMap,WTW,spSz,seMap0,seSel,smoBase)
% sp2graph convert super pixels to curves and graph nodes for GTW
% Input super pixels are not perfect:
% Empty or too small
% Bad corresponding curves
mIhw = major.ihw;
if(isempty(mIhw))
    mIhw = [];
    for i = 1:numel(spLst)
       mIhw = [mIhw;spLst{i}] ;
    end
end
peakTW = major.TW;
tPeak = major.tPeak;
[H0,W0,T0] = size(spMap);
refBase = mean(dFVec(mIhw,:),1);

% make sure the maximum extension.
t00 = min(WTW);
seMap0 = reshape(seMap0,[],T0);
seMap0 = seMap0(mIhw,:);
for t = min(WTW)-1:-1:1
    labels = setdiff(seMap0(:,t),0);
    if(isempty(labels) || ismember(seSel,labels))
        t00 = t;
    else
        break; 
    end
end

ts = extendLeft(imgaussfilt(refBase,2),tPeak,min(peakTW),min(WTW),t00);
ext = 5;
validPart = false(T0,1);
validPart(ts:min(tPeak+ext,T0)) = true;
sz = cellfun(@numel,spLst);

%% make it non-decreasing
smo = 1;
while(1)
    curveSmo = imgaussfilt(refBase,smo);
    left = false(T0,1);
    left(1) = true;
    left(2:end) = curveSmo(2:end)>=curveSmo(1:end-1);
    right = false(T0,1);
    right(end) = true;
    right(1:end-1) = curveSmo(1:end-1)>=curveSmo(2:end);
    maxima = left & right & validPart;
    maxima = find(maxima);
    if(numel(maxima)<=1 || smo>=8)
       break;
    end
    smo = smo + 1;
end
ts = extendLeft(curveSmo,tPeak,min(peakTW),min(WTW),t00);
[~,t00] = min(curveSmo(ts:tPeak)); ts = t00 + ts - 1;
refBase = curveSmo;
peakTW = ts:min(major.tPeak+ext,T0);
refBase(1:ts) = refBase(ts);
refBase(max(peakTW):end) = refBase(max(peakTW));
if(~isempty(maxima))
    [~,tPeak] = max(curveSmo(peakTW));
    tPeak = ts + tPeak - 1;
end
refBase = refBase - refBase(1);
refBase = refBase/refBase(tPeak);
refBase(tPeak:end) = 1;

%% make reference curve non decreasing
refBase = monoIncreasing(refBase,tPeak,ts);

% initialization
nSp = numel(spLst);
tst = zeros(nSp,T0);
ref = zeros(nSp,T0);
orgCurves = zeros(nSp,T0);
stds = zeros(nSp,1);
smoBases = ones(nSp,1)*smoBase;

%% neighbor relation
% graph, at most one pair between two nodes
s = nan(nSp,1);
t = nan(nSp,1);
riseDur = (tPeak-ts+1)*ones(nSp,1);
majorMap = false(H0,W0);
majorMap(mIhw) = true;
inMajor = false(nSp,1);
neiLst = cell(nSp,1);
nPair = 0;
dh = [-1 0 1 -1 0 1 -1 0 1];
dw = [-1 -1 -1 0 0 0 1 1 1];
spMap1 = zeros(H0,W0);
for ii=1:numel(spLst)
    spMap1(spLst{ii}) = ii;
end
nSp = numel(spLst);
for ii=1:nSp
    sp0 = spLst{ii};
    inMajor(ii) = sum(majorMap(sp0)>0);
    [ih0,iw0] = ind2sub([H0,W0],sp0);
    neib0 = [];
    for jj=1:numel(dh)  % find neighbors in eight directions
        ih = max(min(H0,ih0+dh(jj)),1);
        iw = max(min(W0,iw0+dw(jj)),1);
        ihw = sub2ind([H0,W0],ih,iw);
        neib0 = union(neib0,spMap1(ihw));
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

% margin super pixel, add smoothness
for i = 1:nSp
   if(neiLst{i}<3) 
      smoBases(i)  = smoBase*2;
   end
end

%% find seed super pixel
mTW = major.TW;
mtPeak = major.tPeak;
estPeakShift = max(5,round(numel(major.TW)/10));
curves = zeros(nSp,numel(major.TW));
for ii = 1:nSp
    if(inMajor(ii))
        sp0 = spLst{ii};
        curves(ii,:) = mean(dFVec(sp0,major.TW),1);
    end
end
mCurve = mean(dFVec(mIhw,major.TW),1);
cor = corr(mCurve',curves');
[~,seedId] = max(cor);
peakPos = zeros(nSp,1);

visited = false(nSp,1);
checkList = [];
%% seed peak
range = max(1,mtPeak-estPeakShift):min(mtPeak+estPeakShift,T0);
curCurve = imgaussfilt(mean(dFVec(spLst{seedId},:),1),2);
if(curCurve(tPeak)<=curCurve(ts))
    curCurve = mean(dFVec(spLst{seedId},:),1);
end
orgCurves(seedId,:) = curCurve;
curPeak = getPeakPos(curCurve,range,T0,mtPeak,estPeakShift);
visited(seedId) = true;
peakPos(seedId) = curPeak;
ts1 = max(min(peakTW)-ext,1); [~,ts] = min(curCurve(ts1:tPeak)); ts = ts + ts1 - 1;
stds(seedId) =  mean((curCurve(2:end)-curCurve(1:end-1)).^2)/2;
vMin = curCurve(ts);vMax = curCurve(curPeak);
curCurve(1:ts-1) = vMin;
curCurve(curPeak+1:end) = vMax;
curCurve = curCurve - vMin;
curCurve = monoIncreasing(curCurve,curPeak,ts);
tst(seedId,:) = curCurve;
ref(seedId,:) = refBase*(vMax-vMin);
checkList = union(checkList,neiLst{seedId});

while(~isempty(checkList))
   for i = 1:numel(checkList) 
        curLabel = checkList(i);
        sp0 = spLst{curLabel};
        curCurve = imgaussfilt(mean(dFVec(sp0,:),1),2);
        orgCurves(curLabel,:) = curCurve;
        neib0 = neiLst{curLabel};
        neiVisit = neib0(visited(neib0));
        range = max(1,min(peakPos(neiVisit))-estPeakShift):min(max(peakPos(neiVisit))+estPeakShift,T0);
        curPeak = getPeakPos(curCurve,range,T0,round(mean(peakPos(neiVisit))),estPeakShift);
        visited(curLabel) = true;
        peakPos(curLabel) = curPeak;
        ts1 = max(min(peakTW)-ext,1); [~,ts] = min(curCurve(ts1:tPeak)); ts = ts + ts1 - 1;
        stds(curLabel) =  mean((curCurve(2:end)-curCurve(1:end-1)).^2)/2;
        vMin = curCurve(ts);vMax = curCurve(curPeak);
        curCurve(1:ts-1) = vMin;
        curCurve(curPeak+1:end) = vMax;
        curCurve = curCurve - vMin;
        curCurve = monoIncreasing(curCurve,curPeak,ts);
        tst(curLabel,:) = curCurve;
        ref(curLabel,:) = refBase*(vMax-vMin);
        checkList = union(checkList,neiLst{curLabel});
   end
   checkList = checkList(~visited(checkList));
end

idxGood = var(tst,0,2)>1e-10;

end
function [curve,tPeak] = balanceCurve(curve)
    [vP,tPeak] = max(curve);
    [vM,tMin] = min(curve(1:tPeak));
    curve(1:tMin) = vM;
    curve(tPeak:end) = vP;
    curve = curve - vM;
    curve = curve/max(curve);
end
function ts = extendLeft(x0Smo,tPeak,tMin,t00,t01)
%     x0 = x0Smo;
% 	s0 = mean((x0Smo(2:end)-x0Smo(1:end-1)).^2)/2;
    
    [vMin,~] = min(x0Smo(t00:tPeak));
    thr = 0.3*(x0Smo(tPeak)-vMin);
%     tMin = tPeak; vMin = x0Smo(tPeak);
    ts = tMin;
    vMin = x0Smo(tMin);
    
    %% extend backward
    for t = tMin:-1:t01
        if(x0Smo(t)-vMin>=thr)
            break;
        else
           if(x0Smo(t)<vMin) 
              vMin = x0Smo(t) ;
              ts = t;
           end
        end
    end
end
function refBase = monoIncreasing(refBase,tPeak,ts)
    preV = refBase(tPeak);
    preT = tPeak;
    for t = tPeak:-1:ts
       if(refBase(t)<preV) 
           refBase(t:preT) = refBase(t) + (preV-refBase(t))/(preT-t)*[0:preT-t];
           preV = refBase(t);
           preT = t;
       end
    end
end
function maxima = getMaxima(x0Smo,valid,T0)
    curTW = max(min(valid)-2,1):min(T0,max(valid)+2);
    left = false(numel(curTW),1);  left(2:end) = x0Smo(min(curTW)+1:max(curTW))>=x0Smo(min(curTW):max(curTW)-1);
    right = false(numel(curTW),1); right(1:end-1) = x0Smo(min(curTW):max(curTW)-1)>=x0Smo(min(curTW)+1:max(curTW));
    maxima = left & right; maxima = find(maxima); maxima = maxima + min(curTW) - 1;
end
function tPeak = getPeakPos(x0Smo,valid,T0,mtPeak,estPeakShift)
    maxima = getMaxima(x0Smo,valid,T0);
    score = (x0Smo(maxima)').*exp(-(maxima-mtPeak).^2/2/estPeakShift^2/9);
    [~,id] = max(score);
%     [~,id] = max(x0Smo(maxima));
    tPeak = maxima(id);
    if(isempty(maxima))
        tPeak = mtPeak;
    end
end